# Simulation test to test the internal meat of sample_forward.

library(BDAepimodel)
library(mcmc)
library(ggplot2)
library(reshape2)
library(batch)
library(MASS)
library(MCMCpack)


# Initialize epimodel -----------------------------------------------------

set.seed(52786)
require(compiler, quietly = TRUE)

compilePKGS(enable=TRUE)
enableJIT(3)

popsize <- 10

r_meas_process <- function(state, meas_vars, params){
          rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
}

d_meas_process <- function(state, meas_vars, params, log = TRUE) {
          dbinom(x = state[, paste0(meas_vars, "_observed")], size = state[, paste0(meas_vars, "_augmented")], prob = params["rho"], log = log) 
}


epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.25),
                          popsize = popsize,
                          states = c("S", "I", "R"), 
                          params = c(beta = rnorm(1, 0.5, 1e-5), mu = rnorm(1, 1, 1e-5), rho = 0.5, S0 = 0.8, I0 = 0.2, R0 = 0), 
                          rates = c("beta * I", "mu"), 
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                          meas_vars = "I",
                          r_meas_process = r_meas_process,
                          d_meas_process = d_meas_process)


# Simulate the epidemic ---------------------------------------------------

epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = FALSE)


# Set the MCMC settings ---------------------------------------------------

to_estimation_scale <- function(params, epimodel) {
          
          params_est          <- params
          
          params_est["beta"]  <- log(params["beta"] * epimodel$popsize/params["mu"])
          
          params_est["mu"]    <- log(params["mu"])
          
          return(params_est)
}

from_estimation_scale <- function(params_est, epimodel) {
          
          params              <- params_est
          
          params["beta"]      <- exp(params_est["mu"])/epimodel$popsize * exp(params_est["beta"])
          
          params["mu"]        <- exp(params_est["mu"])
          
          return(params)
}

beta_mu_kernel <- function(epimodel) {
          
          # get existing parameter values 
          proposal <- epimodel$params
          
          # get prior likelihood and complete data log likelihood for current parameters
          prior_prob_cur <- c(dnorm(epimodel$params["beta"], mean = log(1), sd = 1.8, log = TRUE),
                              dnorm(epimodel$params["mu"], mean = log(1), sd = 3, log = TRUE)) 
          
          log_likelihood_cur <- epimodel$likelihoods$pop_likelihood_cur + epimodel$likelihoods$obs_likelihood
          
          # transform to estimation scale
          proposal[c("beta", "mu")] <- epimodel$sim_settings$to_estimation_scale(proposal[c("beta", "mu")], epimodel)
          
          # make proposal
          proposal[c("beta", "mu")] <- mvrnorm(n = 1, mu = proposal[c("beta", "mu")], Sigma = epimodel$sim_settings$cov_mtx)
          
          # get prior likelihood for proposed parameters
          prior_prob_new <- c(dnorm(proposal[[1]], mean = log(1), sd = 1.8, log = TRUE),
                              dnorm(proposal[[2]], mean = log(1), sd = 3, log = TRUE))
          
          # transform back to the model scale
          proposal[c("beta","mu")] <- epimodel$sim_settings$from_estimation_scale(proposal[c("beta", "mu")], epimodel)
          
          
          # get population level complete data likelihood for the new parameters
          pop_likelihood_new  <- calc_pop_likelihood(epimodel = epimodel, params = proposal, log = TRUE) 
          log_likelihood_new  <- pop_likelihood_new + epimodel$likelihoods$obs_likelihood
          
          # compute the acceptance probability and accept or reject proposal
          accept_prob <- log_likelihood_new - log_likelihood_cur + 
                    prior_prob_new - prior_prob_cur
          
          if(accept_prob >= 0 || accept_prob >= log(runif(1))) {
                    update_params(epimodel, params = proposal, pop_likelihood = pop_likelihood_new)
                    
          } 
}

rho_kernel <- function(epimodel) {
          
          # Gibbs updates for rho - beta(1, 1) prior
          params_new          <- epimodel$params
          params_new["rho"]   <- rbeta(1, shape1 = 1 + sum(epimodel$obs_mat[,"I_observed"]), shape2 = 1 + sum(epimodel$obs_mat[,"I_augmented"] - epimodel$obs_mat[,"I_observed"]))
          
          obs_likelihood_new <- calc_obs_likelihood(epimodel, params = params_new, log = TRUE)
          
          update_params(epimodel, params = params_new, obs_likelihood = obs_likelihood_new)
          
}

# save new parameters every iteration
# save every tenth configuration matrix
# resample 10 subject-level trajectories in between parameter updates
epimodel$sim_settings <- init_settings(niter = 250000,
                                       save_params_every = 1, 
                                       save_configs_every = 1,
                                       kernel = list(beta_mu_kernel, rho_kernel),
                                       cov_mtx = diag(c(0.005, 0.1), nrow = 2, ncol = 2),
                                       configs_to_redraw = 10,
                                       to_estimation_scale = to_estimation_scale,
                                       from_estimation_scale = from_estimation_scale)

# objects to store results
BDAepimodel_results <- matrix(0, nrow = length(epimodel$obstimes), ncol = epimodel$sim_settings$niter)
Gillespie_results <- matrix(0, nrow = length(epimodel$obstimes), ncol = epimodel$sim_settings$niter)
log_likelihoods <- rep(0, length = epimodel$sim_settings$niter)


# convert epimodel list to environment to enable multiple assignment
.epimodel <- list2env(epimodel, parent = emptyenv(), hash = TRUE)

# move simulation settings to internal objects
niter               <- .epimodel$sim_settings$niter
burnin              <- .epimodel$sim_settings$burnin
save_params_every   <- .epimodel$sim_settings$save_params_every
save_configs_every  <- .epimodel$sim_settings$save_configs_every
kernel              <- .epimodel$sim_settings$kernel
cov_mtx             <- .epimodel$sim_settings$cov_mtx
configs_to_redraw   <- .epimodel$sim_settings$configs_to_redraw
config_replacement  <- ifelse(configs_to_redraw > .epimodel$popsize, TRUE, FALSE)
to_estimation_scale <- .epimodel$sim_settings$to_estimation_scale
from_estimation_scale <- .epimodel$sim_settings$from_estimation_scale
seed                <- .epimodel$sim_settings$seed

# initialize the list of log-likelihoods
.epimodel$likelihoods <- list(pop_likelihood_cur = NULL,
                              pop_likelihood_new = NULL,
                              subj_likelihood_cur = NULL,
                              subj_likelihood_new = NULL,
                              obs_likelihood  = NULL)

# identify whether there is structure in the flow matrix that can be
# leveraged when resampling subject-level trajectories
detect_structure(.epimodel)

# add buffer to the configuration matrix, also adds configurations at
# observation times and instatiates .epimodel$.ind_final_config
.epimodel$config_mat <- expand_config_mat(.epimodel)

# Initialize two lists for the transition probability matrices, one 
# for the matrices and one for products. Also initialize a bookkeeping
# vector indicating which tpms and tpm product matrices need to be 
# recomputed. All tpms and products need to be recomputed following 
# parameter updates, but only the tpms and products for the intervals 
# affected by a change in one of the compartments indicated by 
# index_states need to be recomputed when a new trajectory is redrawn.
.epimodel$.tpms               <- vector("list", length = nrow(.epimodel$config_mat))
.epimodel$.tpm_products       <- vector("list", length = nrow(.epimodel$config_mat))

# initialize the matrix of emission probabilities
.epimodel$nobs                <- length(.epimodel$obstimes)
.epimodel$.emission_mat       <- matrix(0, nrow = .epimodel$num_states, ncol = .epimodel$nobs, dimnames = list(.epimodel$states, .epimodel$obstimes))

# initialize the forward-backward matrices
.epimodel$.fb_mats             <- vector("list", length = .epimodel$nobs - 1)

# initialize the vector of subject-level initial state probabilities
.epimodel$.initdist <- build_initdist(.epimodel)

# get indices for the subject configuration portion of the configuration matrix
.epimodel$.config_inds <- which(grepl(".X", colnames(.epimodel$config_mat)))

# compute the population level likelihood, and the measurement
# process likelihood
.epimodel$likelihoods$pop_likelihood_cur <- calc_pop_likelihood(epimodel = .epimodel, log = TRUE)
.epimodel$likelihoods$obs_likelihood <- calc_obs_likelihood(.epimodel, log = TRUE)

# initialize list for storing results
results <- init_results(.epimodel)

# set seed
results$seed <- seed
set.seed(seed)

# store the initial values in the results matrix
results$params[1, ] <- .epimodel$params
results$log_likelihood[1] <- .epimodel$likelihoods$obs_likelihood + .epimodel$likelihoods$pop_likelihood_cur

results$configs[[1]] <- .epimodel$config_mat[complete.cases(.epimodel$config_mat),]
# initialize the list of log-likelihoods
.epimodel$likelihoods <- list(pop_likelihood_cur = NULL,
                              pop_likelihood_new = NULL,
                              subj_likelihood_cur = NULL,
                              subj_likelihood_new = NULL,
                              obs_likelihood  = NULL)


# identify whether there is structure in the flow matrix that can be
# leveraged when resampling subject-level trajectories
detect_structure(.epimodel)

# add buffer to the configuration matrix, also adds configurations at
# observation times and instatiates .epimodel$.ind_final_config
.epimodel$config_mat <- expand_config_mat(.epimodel)

# Initialize two lists for the transition probability matrices, one 
# for the matrices and one for products. Also initialize a bookkeeping
# vector indicating which tpms and tpm product matrices need to be 
# recomputed. All tpms and products need to be recomputed following 
# parameter updates, but only the tpms and products for the intervals 
# affected by a change in one of the compartments indicated by 
# index_states need to be recomputed when a new trajectory is redrawn.
.epimodel$.tpms               <- vector("list", length = nrow(.epimodel$config_mat))
.epimodel$.tpm_products       <- vector("list", length = nrow(.epimodel$config_mat))

# initialize the matrix of emission probabilities
.epimodel$nobs                <- length(.epimodel$obstimes)
.epimodel$.emission_mat       <- matrix(0, nrow = .epimodel$num_states, ncol = .epimodel$nobs, dimnames = list(.epimodel$states, .epimodel$obstimes))

# initialize the forward-backward matrices
.epimodel$.fb_mats             <- vector("list", length = .epimodel$nobs - 1)

# initialize the vector of subject-level initial state probabilities
.epimodel$.initdist <- build_initdist(.epimodel)

# get indices for the subject configuration portion of the configuration matrix
.epimodel$.config_inds <- which(grepl(".X", colnames(.epimodel$config_mat)))

# compute the population level likelihood, and the measurement
# process likelihood
.epimodel$likelihoods$pop_likelihood_cur <- calc_pop_likelihood(epimodel = .epimodel, log = TRUE)

log_likelihoods <- rep(0, niter)
log_likelihoods[1] <- .epimodel$likelihoods$pop_likelihood_cur 

# initialize the object for drawing the marginal distribution of trajectories directly
epimodel_gillespie <- as.list(.epimodel, all.names = TRUE)

start.time <- Sys.time()


# Start Geweke test -------------------------------------------------------

for(k in 1:niter) {

          if(k%%1000 == 0) print(k)
          
          # choose which subjects should be redrawn
          .subjects <- sample.int(n = .epimodel$popsize, size = configs_to_redraw, replace = config_replacement)
          
          # cycle through subject-level trajectories to be re-drawn
          for(j in 1:configs_to_redraw) {
                   
                    # check to see if any additional irms are needed.
                    # if so, check_irm will instatiate the required
                    # matrices and their eigen decompositions
                    check_irm(.epimodel)
                    
                    # remove trajectory from the counts in config_mat 
                    # and obs_mat, and update the tpm sequences to 
                    # reflect the removal. 
                    remove_trajectory(.epimodel, subject = .subjects[j], save_path = TRUE)
                    
                    # update instatiate missing IRMs, update the tpms 
                    # and tpm products, the emission probability mtx,
                    # and the FB matrices.
                    update_matrices(.epimodel, subject = .subjects[j])
                    
                    # draw a new subject level trajectory
                    draw_trajec(.epimodel, subject = .subjects[j])
                    
          }
          
          # save count of infecteds
          BDAepimodel_results[,k] <- .epimodel$obs_mat[,"I_augmented"]
          
          # simulate a new epidemic via Gillespie and ensure that it consists of data
          Gillespie_results[,k] <- simulate_epimodel(epimodel_gillespie, return_config = FALSE, trim = FALSE)$obs_mat[,"I_augmented"]
          
          if(all(Gillespie_results[,k] == 0)) {
                    while(all(Gillespie_results[,k] == 0)) {
                              Gillespie_results[,k] <- simulate_epimodel(epimodel_gillespie, return_config = FALSE, trim = FALSE)$obs_mat[,"I_augmented"]
                    }
          }
          
          # redraw binomial samples
          .epimodel$obs_mat[,"I_observed"] <- rbinom(n = .epimodel$nobs, .epimodel$obs_mat[,"I_augmented"], .epimodel$params["rho"])
          
          # ensure that not all binomial counts are 0, i.e. that there is data
          if(all(.epimodel$obs_mat[,"I_observed"]==0)) {
                    while(all(.epimodel$obs_mat[,"I_observed"]==0)) {
                              .epimodel$obs_mat[,"I_observed"] <- rbinom(n = .epimodel$nobs, .epimodel$obs_mat[,"I_augmented"], .epimodel$params["rho"])
                    }  
          }
          
          # save log_likelihood
          log_likelihoods[k] <- .epimodel$likelihoods$pop_likelihood_cur 
          
}
end.time <- Sys.time()


# Save results ------------------------------------------------------------

results <- list(BDAepimodel_results, Gillespie_results, time = difftime(end.time, start.time, units = "mins"))
save(results, file = "BDAepimodel_geweke.Rdata")


# Generate Plots ----------------------------------------------------------

pdf(file = "BDAepimodel_Geweke_plots.pdf")

results_comp <- data.frame(method = rep(c("BDAepimodel", "Gillespie Counts"), each = .epimodel$nobs), time = rep(.epimodel$obstimes, 2), count = 0, mcse = 0)
results_comp[results_comp$method == "BDAepimodel","count"] <- rowMeans(BDAepimodel_results)
results_comp[results_comp$method == "Gillespie Counts","count"] <- rowMeans(Gillespie_results)
BDAepimodel_vars <- sapply(apply(BDAepimodel_results, 1, initseq), "[[", "var.pos")
results_comp[results_comp$method == "BDAepimodel","mcse"] <- sqrt(BDAepimodel_vars)/sqrt(niter)
results_comp[results_comp$method == "Gillespie","mcse"] <- apply(Gillespie_results, 1, sd)/sqrt(niter)

ggplot(results_comp, aes(x = time, y = count, colour = method)) + geom_ribbon(aes(ymin = count - 1.96 * mcse, ymax = count + 1.96 * mcse, fill = method), alpha = 0.2) + geom_line() + theme_bw() + labs(title = "BDAepimodel vs. Gillespie")

dev.off()


save(results, file = "BDAepimodel_geweke.Rdata")
