# this file contains pieces of code that are used in developing the package but
# are not necessary for the package to work.


# Initialize epimodel -----------------------------------------------------

set.seed(52787)
require(BDAepimodel)

popsize <- 10

r_meas_process <- function(state, meas_vars, params){
          rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
}

d_meas_process <- function(state, meas_vars, params, log = TRUE) {
          dbinom(x = state[, paste0(meas_vars, "_observed")], size = state[, paste0(meas_vars, "_augmented")], prob = params["rho"], log = log) 
}



epimodel <- init_epimodel(obstimes = seq(0, 10, by = 1),
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
epimodel$sim_settings <- init_settings(niter = 100,
                                       save_params_every = 1, 
                                       save_configs_every = 1,
                                       kernel = list(beta_mu_kernel, rho_kernel),
                                       cov_mtx = diag(c(0.005, 0.1), nrow = 2, ncol = 2),
                                       configs_to_redraw = 200,
                                       to_estimation_scale = to_estimation_scale,
                                       from_estimation_scale = from_estimation_scale)


# fit the model -----------------------------------------------------------

Rprof(prof <- tempfile())
epimodel <- fit_epimodel(epimodel, monitor = FALSE) 
Rprof()
summaryRprof(prof)

require(proftools)
plotProfileCallGraph(readProfileData(prof), score = "self")
