# this file contains pieces of code that are used in developing the package but
# are not necessary for the package to work.
popsize <- 200

initialization_fcn <- function(){
          init_state <- as.numeric(rmultinom(1, popsize, c(0.95, 0.05, 0)))
          names(init_state ) <- c("S", "I", "R")
          return(init_state)
}

r_meas_process <- function(state, meas_vars, params){
          rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
}

d_meas_process <- function(state, meas_vars, params, log = TRUE) {
          dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_truth", sep = "")], prob = params["rho"], log = log) 
}

# evaluates initial distribution for a single subject
d_initdist <- function(state, params, log = TRUE) {
          if(log == TRUE) {
                    (state == 2)* log(params["p0"]) + (state == 1) * log(1-params["p0"])
          } else {
                    (params["p0"] ^ (state == 2)) * ((1-params["p0"])^(state == 1))
          }
}

# subject level simulation of initial state at time t0
r_initdist <- function(params) {
          sample.int(3, 1, prob = c(1-params["p0"], params["p0"], 0))
          }

# R0 = 4, mu = 1, rho = 0.5, p0 = 0.05
epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.5),
                          popsize = popsize,
                          states = c("S", "I", "R"), 
                          params = c(beta = 0.02, mu = 1, rho = 0.5, p0 = 0.05), 
                          rates = c("beta * I", "mu"), 
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                          meas_vars = "I",
                          r_meas_process = r_meas_process,
                          d_meas_process = d_meas_process,
                          d_initdist = d_initdist,
                          r_initdist = r_initdist)

epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = FALSE)


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
                                                         [,"I_observed"]), shape2 = 1 + sum(epimodel$obs_mat[,"I_truth"] - 
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
epimodel$sim_settings <- init_settings(niter = 1000,
                burnin <- 200,
                params_every <- 1, 
                configs_every <- 10,
                kernel = kernel,
                cov_mtx = cov_mtx,
                configs_to_redraw = 20,
                to_estimation_scale = to_estimation_scale,
                from_estimation_scale = from_estimation_scale)