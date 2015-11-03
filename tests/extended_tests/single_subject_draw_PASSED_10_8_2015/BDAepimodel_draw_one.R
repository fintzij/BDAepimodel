# Simulation test to test the internal meat of sample_forward.

library(BDAepimodel)
library(mcmc)
library(GillespieSSA)
library(doParallel)
library(doMC)
library(itertools)
library(foreach)
library(iterators)
library(ggplot2)
library(reshape2)
library(batch)

parseCommandArgs()

# Set up simulation objects ----------------------------------------------------------------

# this file contains pieces of code that are used in developing the package but
# are not necessary for the package to work.
set.seed(52787)
require(BDAepimodel)

popsize <- 10

r_meas_process <- function(state, meas_vars, params){
          rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
}

d_meas_process <- function(state, meas_vars, params, log = TRUE) {
          dbinom(x = state[, paste0(meas_vars, "_observed")], size = state[, paste0(meas_vars, "_augmented")], prob = params["rho"], log = log) 
}

epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.05),
                          popsize = popsize,
                          states = c("S", "I", "R"), 
                          params = c(beta = rnorm(1, 0.5, 1e-5), mu = rnorm(1, 1, 1e-5), rho = 0.5, S0 = 0.8, I0 = 0.2, R0 = 0), 
                          rates = c("beta * I", "mu"), 
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                          meas_vars = "I",
                          r_meas_process = r_meas_process,
                          d_meas_process = d_meas_process)

epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = FALSE)
epimodel_gillespie <- epimodel



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
          proposal[c("beta", "mu")] <- MASS::mvrnorm(n = 1, mu = proposal[c("beta", "mu")], Sigma = epimodel$sim_settings$cov_mtx)
          
          # get prior likelihood for proposed parameters
          prior_prob_new <- c(dnorm(proposal[[1]], mean = log(1), sd = 1.8, log = TRUE),
                              dnorm(proposal[[2]], mean = log(1), sd = 3, log = TRUE))
          
          # transform back to the model scale
          proposal[c("beta","mu")] <- epimodel$sim_settings$from_estimation_scale(proposal[c("beta", "mu")], epimodel)
          
          
          # update array of rate matrices
          epimodel <- build_new_irms(epimodel, proposal)
          
          # get population level complete data likelihood for the new parameters
          pop_likelihood_new  <- calc_pop_likelihood(epimodel, log = TRUE)
          obs_likelihood_new  <- calc_obs_likelihood(epimodel = epimodel, params = proposal, log = TRUE)
          log_likelihood_new  <- pop_likelihood_new + obs_likelihood_new
          
          # compute the acceptance probability and accept or reject proposal
          accept_prob <- log_likelihood_new - log_likelihood_cur + sum(prior_prob_new) - sum(prior_prob_cur)
          
          if(accept_prob >= 0 || accept_prob >= log(runif(1))) {
                    
                    # update parameters, likelihood objects, and eigen decompositions
                    epimodel <- update_params(epimodel, params = proposal, pop_likelihood = pop_likelihood_new, obs_likelihood = obs_likelihood_new)
                    
                    # update the eigen decompositions
                    buildEigenArray(eigenvals = epimodel$eigen_values, eigenvecs = epimodel$eigen_vectors, inversevecs = epimodel$inv_eigen_vectors, irm_array = epimodel$irm)
                    
          } else {
                    # restore the previous array of rate matrices
                    epimodel$irm <- epimodel$.irm_cur
          }
          
          return(epimodel)
}

rho_kernel <- function(epimodel) {
          
          # Gibbs updates for rho - beta(1, 1) prior
          params_new          <- epimodel$params
          params_new["rho"]   <- rbeta(1, shape1 = 1 + sum(epimodel$obs_mat[,"I_observed"]), shape2 = 1 + sum(epimodel$obs_mat[,"I_augmented"] - epimodel$obs_mat[,"I_observed"]))
          
          obs_likelihood_new <- calc_obs_likelihood(epimodel, params = params_new, log = TRUE)
          
          epimodel <- update_params(epimodel, params = params_new, obs_likelihood = obs_likelihood_new)
          
          return(epimodel)
}

kernel <- list(beta_mu_kernel, rho_kernel)


# save new parameters every iteration
# save every tenth configuration matrix
# resample 10 subject-level trajectories in between parameter updates
epimodel <- init_settings(epimodel,
                          niter = 250000,
                          save_params_every = 1, 
                          save_configs_every = 1,
                          kernel = kernel,
                          cov_mtx = diag(c(0.02, 0.02, 0.02), nrow = 3, ncol = 3),
                          configs_to_redraw = 10,
                          to_estimation_scale = to_estimation_scale,
                          from_estimation_scale = from_estimation_scale)



# objects to store results
BDAepimodel_results <- matrix(0, nrow = length(epimodel$obstimes), ncol = epimodel$sim_settings$niter)
Gillespie_results <- matrix(0, nrow = length(epimodel$obstimes), ncol = epimodel$sim_settings$niter)


epimodel <- prepare_epimodel(epimodel)

# move some of the simulation settings to internal objects
niter                 <- epimodel$sim_settings$niter
save_params_every     <- epimodel$sim_settings$save_params_every
save_configs_every    <- epimodel$sim_settings$save_configs_every
kernel                <- epimodel$sim_settings$kernel
configs_to_redraw     <- epimodel$sim_settings$configs_to_redraw
config_replacement    <- epimodel$sim_settings$config_replacement

# initialize list for storing results
results <- init_results(epimodel)

# set seed
results$seed        <- epimodel$sim_settings$seed
set.seed(results$seed)

# store the initial values in the results matrix
results$params[1, ] <- epimodel$params
results$log_likelihood[1] <- epimodel$likelihoods$obs_likelihood + epimodel$likelihoods$pop_likelihood_cur

results$configs[[1]] <- epimodel$pop_mat[1:epimodel$ind_final_config,]

# initialize the vector for path acceptances
epimodel$path_accept_vec <- rep(0, configs_to_redraw)

# get start time
start.time = Sys.time()

for(k in 1:niter) {
          
          if(k%%500 == 0) print(k)
          
         # resimulate the epidemic
          .epimodel <- simulate_epimodel(epimodel, trim = FALSE, lump = TRUE)
          if(.epimodel$pop_mat[1,"I"] < 1) {
                   while(.epimodel$pop_mat[1,"I"] < 1) {
                             .epimodel <- simulate_epimodel(epimodel, trim = FALSE, lump = TRUE)
                   } 
          }
          
          .epimodel<- prepare_epimodel(.epimodel)
          
          if(all(.epimodel$obs_mat[,"I_observed"]==0)) {
                  while(all(.epimodel$obs_mat[,"I_observed"]==0)) {
                            .epimodel$obs_mat[,"I_observed"] <- rbinom(n = .epimodel$nobs, .epimodel$obs_mat[,"I_augmented"], .epimodel$params["rho"])
                  }  
          }
          
          # choose which subjects should be redrawn
          subjects <- sample.int(10, 1)
          
          # save the original trajectory to the gillespie_results object
          retrieveSubjPath(.epimodel$subj_path, subjects, .epimodel$pop_mat, .epimodel$init_config, .epimodel$ind_final_config, .epimodel$flow_inds)
          Gillespie_results[,k] <- ifelse(.epimodel$subj_path[.epimodel$obs_time_inds] == 2, 1, 0)
          
          # remove trajectory from the counts in config_mat 
          # and obs_mat, and update the tpm sequences to 
          # reflect the removal. 
          .epimodel <- remove_trajectory(.epimodel, subject = subjects, save_path = TRUE)
          
          
          # TPM sequence
          tpmSeqs(tpms = .epimodel$tpms, pop_mat = .epimodel$pop_mat, eigen_vals = .epimodel$eigen_values, eigen_vecs = .epimodel$eigen_vectors, inverse_vecs = .epimodel$inv_eigen_vectors, irm_keys = .epimodel$keys)
          
          # TPM product subsequences
          tpmProdSeqs(tpm_prods = .epimodel$tpm_products, tpms = .epimodel$tpms, obs_time_inds = .epimodel$obs_time_inds)
          
          # Emission matrix
          .epimodel$emission_mat <- build_emission_mat(emission_mat = .epimodel$emission_mat, epimodel = .epimodel)
          
          # FB matrices
          buildFBMats(.epimodel$fb_mats, .epimodel$tpm_products, .epimodel$emission_mat, .epimodel$initdist, .epimodel$obs_time_inds)                        
          
          # sample status at observation times
          .epimodel$subj_path <- sample_at_obs_times(path = .epimodel$subj_path, epimodel = .epimodel)
          
          # sample status at event times
          .epimodel$subj_path <- sample_at_event_times(path = .epimodel$subj_path, epimodel = .epimodel)
          
          # sample paths in inter-event intervals
          path <- sample_path(.epimodel, subject = subjects)
          .epimodel$n_jumps <- ifelse(is.null(path), 0, nrow(path))
          
          # insert the path
          if(!is.null(path)) {
                    
                    # add rows to the population level bookkeeping matrix if necessary
                    if(.epimodel$ind_final_config + .epimodel$n_jumps > nrow(.epimodel$pop_mat)) {
                              .epimodel <- expand_pop_mat(.epimodel, buffer_size = .epimodel$n_jumps)
                    }
                    
                    # insert the path
                    insertPath(path, subjects, .epimodel$pop_mat, .epimodel$subj_path, .epimodel$ind_final_config)
          } 
          
          # insert the proposed trajectory
          .epimodel <- insert_trajectory(epimodel = .epimodel, subject = subjects, reinsertion = FALSE)
          # accept or reject the proposed path
          .epimodel <- MH_accept_reject(epimodel = .epimodel, subject = subjects, iter = k)
          
          BDAepimodel_results[,k] <- ifelse(.epimodel$subj_path[.epimodel$obs_time_inds] == 2, 1, 0)
          
}
end.time <- Sys.time()

results <- list(BDAepimodel_results, Gillespie_results, time = difftime(end.time, start.time, units = "mins"))

save(results, file = "BDAepimodel_draw_one.Rdata")


pdf(file = "BDAepimodel_draw_one_plots.pdf")

results_comp <- data.frame(method = rep(c("BDAepimodel", "Gillespie Counts"), each = .epimodel$nobs), time = rep(.epimodel$obstimes, 2), count = 0, mcse = 0)
results_comp[results_comp$method == "BDAepimodel","count"] <- rowMeans(BDAepimodel_results)
results_comp[results_comp$method == "Gillespie Counts","count"] <- rowMeans(Gillespie_results)
BDAepimodel_vars <- apply(BDAepimodel_results, 1, sd)
results_comp[results_comp$method == "BDAepimodel","mcse"] <- sqrt(BDAepimodel_vars)/sqrt(niter)
results_comp[results_comp$method == "Gillespie","mcse"] <- apply(Gillespie_results, 1, sd)/sqrt(niter)

ggplot(results_comp, aes(x = time, y = count, colour = method)) + geom_ribbon(aes(ymin = count - 1.96 * mcse, ymax = count + 1.96 * mcse, fill = method), alpha = 0.2) + geom_line() + theme_bw() + labs(title = "BDAepimodel vs. Gillespie")

dev.off()


save(results, file = "fixedpar_convcomp.Rdata")