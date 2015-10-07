# Simulation test to test the internal meat of sample_forward.

library(BDAepimodel)
library(mcmc)
library(ggplot2)
library(reshape2)
library(batch)

parseCommandArgs()

# Set up simulation objects ----------------------------------------------------------------

# this file contains pieces of code that are used in developing the package but
# are not necessary for the package to work.
set.seed(52787)

popsize <- 20

r_meas_process <- function(state, meas_vars, params){
          rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
}

d_meas_process <- function(state, meas_vars, params, log = TRUE) {
          dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
}

epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.1),
                          popsize = popsize,
                          states = c("S", "I", "R"), 
                          params = c(beta = rnorm(1, 0.25, 1e-6), mu = rnorm(1, 1, 1e-6), rho = 0.5, S0 = 0.8, I0 = 0.2, R0 = 0), 
                          rates = c("beta * I", "mu"), 
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                          meas_vars = "I",
                          r_meas_process = r_meas_process,
                          d_meas_process = d_meas_process)

epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = FALSE)
epimodel_gillespie <- epimodel


to_estimation_scale <- list(function(params, popsize) {log(params["beta"] * popsize/params["mu"])}, 
                            function(params) {log(epimodel$params["mu"])})

from_estimation_scale <- list(function(params_est, popsize) {exp(params_est["mu"])/popsize * exp(params_est["beta"])},
                              function(params_est) {exp(params_est["mu"])})


beta_mu_kernel <- function(epimodel, cov_mtx, log_likelihood_cur, 
                           prior_prob_cur, dprior, to_estimation_scale, from_estimation_scale) {
          
          # four vectors - existing parameters on natural and estimation scale, 
          # and new parameters on natural and estimation scale. 
          params_est <- epimodel$params
          params_est[c("beta","mu")] <- to_estimation_scale(params_est[c("beta",
                                                                         "mu")], epimodel$popsize) 
          
          proposal_est <- proposal <- params_est
          
          # make proposal
          proposal_est[c("beta", "mu")] <- mvrnorm(n = 1, mu = params[c("beta", 
                                                                        "mu")], Sigma = cov_mtx)
          
          proposal[c("beta","mu")]<-from_estimation_scale(proposal_est[c("beta",
                                                                         "mu")], epimodel$popsize)
          
          # update prior probability on estimation scale
          prior_prob_new <- update_prior_probs(params_new = proposal_est, 
                                               params_cur = params_est, dprior = dprior, log = TRUE)
          log_likelihood_new <- calc_log_likelihood(epimodel = epimodel, params 
                                                    = proposal)
          
          # compute the acceptance probability and accept or reject proposal
          accept_prob <- log_likelihood_new - log_likelihood_cur + 
                    prior_prob_new - prior_prob_cur
          
          if(accept_prob >= 0 || accept_prob >= log(runif(1))) {
                    params_new <- proposal
          } else {
                    params_new <- epimodel$params
          }
          
          return(params_new)        
}

rho_p0_kernel <- function(epimodel) {
          
          # Gibbs updates for rho and p0 - beta(1, 1) prior for each
          params_new <- epimodel$params
          params_new["rho"] <- rbeta(1, shape1 = 1 + sum(epimodel$obs_mat
                                                         [,"I_observed"]), shape2 = 1 + sum(epimodel$obs_mat[,"I_augmented"] - 
                                                                                                      epimodel$obs_mat[,"I_observed"]))
          params_new["p0"] <- rbeta(1, shape1 = 1 + epimodel$config_mat[1, "I"],
                                    shape2 = 1 + epimodel$popsize - epimodel$config_mat[1, "I"])
          
          return(params_new)        
}

kernel <- list(beta_mu_kernel, rho_p0_kernel)
cov_mtx <- diag(c(0.005, 0.1), nrow = 2, ncol = 2)
config_resample_prop <- 10

# save new parameters every iteration
# save every tenth configuration matrix
# resample 10 subject-level trajectories in between parameter updates
epimodel$sim_settings <- init_settings(epimodel = epimodel,
                                       niter = 100000,
                                       burnin = 0,
                                       params_every = 1, 
                                       configs_every = 1, 
                                       kernel = kernel,
                                       cov_mtx = cov_mtx,
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
.niter              <- .epimodel$sim_settings$niter
.burnin             <- .epimodel$sim_settings$burnin
.params_every       <- .epimodel$sim_settings$params_every
.configs_every      <- .epimodel$sim_settings$configs_every
.kernel             <- .epimodel$sim_settings$kernel
.cov_mtx            <- .epimodel$sim_settings$cov_mtx
.configs_to_redraw  <- .epimodel$sim_settings$configs_to_redraw
.config_replacement <- ifelse(.configs_to_redraw > .epimodel$popsize, TRUE, FALSE)
.parallelize        <- .epimodel$sim_settings$parallelize
.to_estimation_scale <- .epimodel$sim_settings$to_estimation_scale
.from_estimation_scale <- .epimodel$sim_settings$from_estimation_scale

# initialize the list of log-likelihoods
epimodel$likelihoods <- list(pop_likelihood_cur = NULL,
                             pop_likelihood_new = NULL,
                             subj_likelihood_cur = NULL,
                             subj_likelihood_new = NULL,
                             obs_likelihood  = NULL)

# identify the unmeasured compartments, constants for number of
# measured and unmeasured states
.epimodel$unmeasured_vars     <- setdiff(.epimodel$states, .epimodel$meas_vars)
.epimodel$num_unmeasured      <- length(.epimodel$unmeasured_vars)
.epimodel$num_measured        <- length(.epimodel$meas_vars)
.epimodel$num_states          <- length(.epimodel$states)

# set the indexing for the observed and measured variables in the
# observation matrix
.epimodel$meas_vars_aug       <- paste0(.epimodel$meas_vars, "_augmented")
.epimodel$meas_vars_obs       <- paste0(.epimodel$meas_vars, "_observed")

# identify whether there is structure in the flow matrix that can be
# leveraged when resampling subject-level trajectories
detect_structure(.epimodel)

# initialize list for storing results
.results <- init_results(.epimodel)

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

build_irm(.epimodel)
.epimodel$likelihoods$pop_likelihood_cur <- calc_pop_likelihood(epimodel = .epimodel, log = TRUE)

start.time <- Sys.time()
for(k in 1:.niter) {
          
          if(k%%500 == 0) print(k)
          
          # re-compute the array of rate matrices - only needs to be
          # recomputed every time parameters are updated.
          
          # re-compute the population level likelihood, measurement
          # process likelihood, and complete data log-likelihood
#           .epimodel$likelihoods$pop_likelihood_cur <- calc_pop_likelihood(epimodel = .epimodel, log = TRUE)
          # .epimodel$likelihoods$obs_likelihood_cur <- calc_obs_likelihood(epimodel = .epimodel, log = TRUE)
          
          # choose which subjects should be redrawn
          .subjects <- sample.int(n = .epimodel$popsize, size = .configs_to_redraw, replace = .config_replacement)

          # cycle through subject-level trajectories to be re-drawn
          for(j in 1:.configs_to_redraw) {
                   
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
          
          BDAepimodel_results[,k] <- .epimodel$obs_mat[,"I_augmented"]
          
          Gillespie_results[,k] <- simulate_epimodel(epimodel_gillespie, return_config = FALSE, trim = FALSE)$obs_mat[,"I_augmented"]
           
          .epimodel$obs_mat[,"I_observed"] <- rbinom(n = .epimodel$nobs, .epimodel$obs_mat[,"I_augmented"], .epimodel$params["rho"])
          
          if(all(.epimodel$obs_mat[,"I_observed"]==0)) {
                    while(all(.epimodel$obs_mat[,"I_observed"]==0)) {
                              .epimodel$obs_mat[,"I_observed"] <- rbinom(n = .epimodel$nobs, .epimodel$obs_mat[,"I_augmented"], .epimodel$params["rho"])
                    }  
          }
     
          log_likelihoods[k] <- .epimodel$likelihoods$pop_likelihood_cur     
}
end.time <- Sys.time()

results <- list(BDAepimodel_results, Gillespie_results, time = difftime(end.time, start.time, units = "mins"))
save(results, file = "BDAepimodel_geweke.Rdata")


pdf(file = "BDAepimodel_Geweke_plots.pdf")

results_comp <- data.frame(method = rep(c("BDAepimodel", "Gillespie Counts"), each = .epimodel$nobs), time = rep(.epimodel$obstimes, 2), count = 0, mcse = 0)
results_comp[results_comp$method == "BDAepimodel","count"] <- rowMeans(BDAepimodel_results)
results_comp[results_comp$method == "Gillespie Counts","count"] <- rowMeans(Gillespie_results)
BDAepimodel_vars <- sapply(apply(BDAepimodel_results, 1, initseq), "[[", "var.pos")
results_comp[results_comp$method == "BDAepimodel","mcse"] <- sqrt(BDAepimodel_vars)/sqrt(.niter)
results_comp[results_comp$method == "Gillespie","mcse"] <- apply(Gillespie_results, 1, sd)/sqrt(.niter)

ggplot(results_comp, aes(x = time, y = count, colour = method)) + geom_ribbon(aes(ymin = count - 1.96 * mcse, ymax = count + 1.96 * mcse, fill = method), alpha = 0.2) + geom_line() + theme_bw() + labs(title = "BDAepimodel vs. Gillespie")

dev.off()


save(results, file = "BDAepimodel_geweke.Rdata")
