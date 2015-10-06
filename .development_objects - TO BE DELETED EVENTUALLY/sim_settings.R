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



# R0 = 4, mu = 1, rho = 0.5, p0 = 0.05
epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.5),
                          popsize = popsize,
                          states = c("S", "I", "R"), 
                          params = c(beta = rnorm(1, 0.5, 1e-6), mu = rnorm(1, 1, 1e-6), rho = 0.5, S0 = 0.7, I0 = 0.2, R0 = 0.1), 
                          rates = c("beta * I", "mu"), 
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                          meas_vars = "I",
                          r_meas_process = r_meas_process,
                          d_meas_process = d_meas_process)

epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = FALSE)


to_estimation_scale <- list(function(params, popsize) {log(params["beta"] * popsize/params["mu"])}, 
                            function(params) {log(epimodel$params["mu"])})

from_estimation_scale <- list(function(params_est, popsize) {exp(params_est["mu"])/popsize * exp(params_est["beta"])},
                              function(params_est) {exp(params_est["mu"])})


beta_mu_kernel <- function(epimodel) {
          
          # four vectors - existing parameters on natural and estimation scale, 
          # and new parameters on natural and estimation scale. 
          params_est <- epimodel$params
          params_est[c("beta","mu")] <- to_estimation_scale(params_est[c("beta", "mu")], epimodel$popsize) 
          
          proposal_est <- proposal <- params_est
          
          # make proposal
          proposal_est[c("beta", "mu")] <- mvrnorm(n = 1, mu = params[c("beta", "mu")], Sigma = cov_mtx)
          
          proposal[c("beta","mu")]<-from_estimation_scale(proposal_est[c("beta", "mu")], epimodel$popsize)
          
          # update prior probability on estimation scale
          prior_prob_new      <- update_prior_probs(epimodel, params_new = proposal_est, dprior = dprior, log = TRUE)
          log_likelihood_new  <- calc_log_likelihood(epimodel = epimodel, params = proposal)
          
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

rho_kernel <- function(epimodel) {
          
          # Gibbs updates for rho - beta(1, 1) prior
          params_new          <- epimodel$params
          params_new["rho"]   <- rbeta(1, shape1 = 1 + sum(epimodel$obs_mat[,"I_observed"]), shape2 = 1 + sum(epimodel$obs_mat[,"I_augmented"] - epimodel$obs_mat[,"I_observed"]))
          
          return(params_new)        
}

# save new parameters every iteration
# save every tenth configuration matrix
# resample 10 subject-level trajectories in between parameter updates
epimodel$sim_settings <- init_settings(niter = 1000,
                burnin <- 200,
                params_every <- 1, 
                configs_every <- 10,
                kernel = list(beta_mu_kernel, rho_kernel),
                cov_mtx = diag(c(0.005, 0.1), nrow = 2, ncol = 2),
                configs_to_redraw = 5,
                to_estimation_scale = to_estimation_scale,
                from_estimation_scale = from_estimation_scale)