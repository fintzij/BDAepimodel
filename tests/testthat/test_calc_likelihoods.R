library(BDAepimodel)

context("Calculation of the population-level log-likelihood and measurement processes log-likelihood") 

test_that("The log-likelihoods are computed correctly", {
          
          set.seed(52787)
          
          popsize <- 3
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
          }
          
          # R0 = 4, mu = 1, rho = 0.5, p0 = 0.05
          epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.5),
                                    popsize = popsize,
                                    states = c("S", "I", "R"), 
                                    params =  c(beta = rnorm(1, 1, 1e-6), mu = rnorm(1, 1, 1e-6), rho = 0.5,  S0 = 0.5, I0 = 0.5, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          epimodel <- prepare_epimodel(epimodel)
          
          pop_mat <- epimodel$pop_mat[c(1, which(epimodel$pop_mat[,"ID"] != 0), epimodel$ind_final_config),]
          
          pop_likelihood <- 2 * log(epimodel$params["I0"]) + log( 1- epimodel$params["I0"]) + 
                    log(2 * epimodel$params["beta"]) - (2 + 2)*(pop_mat[2, "time"] - pop_mat[1, "time"]) +
                    log(1) - (3)*(pop_mat[3, "time"] - pop_mat[2, "time"]) +
                    log(1) - (2)*(pop_mat[4, "time"] - pop_mat[3, "time"]) +
                    log(1) - (1)*(pop_mat[5, "time"] - pop_mat[4, "time"]) 
          
          obs_likelihood <-sum(dbinom(epimodel$obs_mat[,"I_observed"], epimodel$obs_mat[,"I_augmented"], prob = epimodel$params["rho"], log= TRUE))
          
          subj1_likelihood <- log(1 - epimodel$params["I0"]) + log(2*epimodel$params["beta"]) - 2*epimodel$params["beta"] * (pop_mat[2,"time"] - pop_mat[1,"time"]) + log(epimodel$params["mu"])  - epimodel$params["mu"] *(pop_mat[3,"time"] - pop_mat[2,"time"])
          
          expect_equal(signif(calc_pop_likelihood(epimodel = epimodel, log = TRUE), digits = 6), as.numeric(signif(pop_likelihood, digits = 6)))
          expect_equal(as.numeric(calc_obs_likelihood(epimodel = epimodel, log = TRUE)), as.numeric(obs_likelihood))
          expect_equal(as.numeric(subjectLikelihood(1, epimodel$pop_mat, epimodel$config_mat, epimodel$irm, epimodel$initdist, epimodel$keys, c(1, c(1:epimodel$ind_final_config)[-epimodel$obs_time_inds], epimodel$ind_final_config), TRUE)), as.numeric(subj1_likelihood))
}) 


test_that("population-level likelihoods are correctly computed within each iteration", {
          
          set.seed(52787)
          require(BDAepimodel)
          
          popsize <- 5 
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste0(meas_vars, "_observed")], size = state[, paste0(meas_vars, "_augmented")], prob = params["rho"], log = log) 
          }
          
          
          
          epimodel <- init_epimodel(obstimes = seq(0, 10, by = 1),
                                    popsize = popsize,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = rnorm(1, 1, 1e-5), mu = rnorm(1, 1, 1e-5), rho = 0.5, S0 = 0.8, I0 = 0.2, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process)
          
          
          # Simulate the epidemic ---------------------------------------------------
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          
          # Set the MCMC settings ---------------------------------------------------
          
          # Kernel - MH - Beta, Mu, Rho ---------------------------------------------
          
          to_estimation_scale <- function(params, epimodel) {
                    
                    params_est          <- params
                    
                    params_est["beta"]  <- log(params["beta"] * epimodel$popsize/params["mu"])
                    
                    params_est["mu"]    <- log(params["mu"])
                    
                    params_est["rho"]   <- logit(params["rho"])
                    
                    return(params_est)
          }
          
          from_estimation_scale <- function(params_est, epimodel) {
                    
                    params              <- params_est
                    
                    params["beta"]      <- exp(params_est["mu"])/epimodel$popsize * exp(params_est["beta"])
                    
                    params["mu"]        <- exp(params_est["mu"])
                    
                    params["rho"]       <- expit(params["rho"])
                    
                    return(params)
          }
          
          
          beta_mu_rho_kernel <- function(epimodel) {
                    
                    # get existing parameter values 
                    proposal <- epimodel$params
                    
                    # get prior likelihood and complete data log likelihood for current parameters
                    prior_prob_cur <- c(dnorm(epimodel$params["beta"], mean = log(1), sd = 1.8, log = TRUE),
                                        dnorm(epimodel$params["mu"], mean = log(1), sd = 3, log = TRUE),
                                        dnorm(epimodel$params["rho"], mean = log(0.4), sd = 0.6, log = TRUE)) 
                    
                    log_likelihood_cur <- epimodel$likelihoods$pop_likelihood_cur + epimodel$likelihoods$obs_likelihood
                    
                    # transform to estimation scale
                    proposal[c("beta", "mu", "rho")] <- epimodel$sim_settings$to_estimation_scale(proposal[c("beta", "mu", "rho")], epimodel)
                    
                    # make proposal
                    proposal[c("beta", "mu", "rho")] <- MASS::mvrnorm(n = 1, mu = proposal[c("beta", "mu", "rho")], Sigma = epimodel$sim_settings$cov_mtx)
                    
                    # get prior likelihood for proposed parameters
                    prior_prob_new <- c(dnorm(proposal[[1]], mean = log(1), sd = 1.8, log = TRUE),
                                        dnorm(proposal[[2]], mean = log(1), sd = 3, log = TRUE),
                                        dnorm(proposal[[3]], mean = log(0.4), sd = 0.6, log = TRUE))
                    
                    # transform back to the model scale
                    proposal[c("beta","mu", "rho")] <- epimodel$sim_settings$from_estimation_scale(proposal[c("beta", "mu", "rho")], epimodel)
                    
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
                              
                    } else {
                              # restore the previous array of rate matrices
                              epimodel$irm <- epimodel$.irm_cur
                    }
                    
                    return(epimodel)
          }
          epimodel <- init_settings(epimodel,
                                    niter = 4,
                                    save_params_every = 1, 
                                    save_configs_every = 5,
                                    kernel = list(beta_mu_rho_kernel),
                                    cov_mtx = diag(c(0.02, 0.02, 0.02), nrow = 3, ncol = 3),
                                    configs_to_redraw = 20,
                                    to_estimation_scale = to_estimation_scale,
                                    from_estimation_scale = from_estimation_scale)
          
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
          
          results$configs[[1]] <- epimodel$pop_mat[complete.cases(epimodel$pop_mat),]
          
          # initialize the vector for path acceptances
          epimodel$path_accept_vec <- rep(0, configs_to_redraw)
          
          # get start time
          start.time = Sys.time()
          
          # generate niter parameter samples
          for(k in 2:niter) {
                    
                    # re-compute the arrays of rate matrices and eigen
                    # decompositions - only needs to be recomputed every time
                    # parameters are updated.
                    
                    # choose which subjects should be redrawn
                    subjects <- sample.int(n = epimodel$popsize, size = configs_to_redraw, replace = config_replacement)
                    
                    # cycle through subject-level trajectories to be re-drawn
                    for(j in 1:configs_to_redraw) {
                              
                              # remove trajectory from the counts in config_mat
                              # and obs_mat, and update the tpm sequences to
                              # reflect the removal.
                              epimodel <- remove_trajectory(epimodel, subject = subjects[j], save_path = TRUE)
                              
                              # update instatiate missing IRMs, update the tpms
                              # and tpm products, the emission probability mtx,
                              # and the FB matrices.
                              
                              # TPM sequence
                              tpmSeqs(tpms = epimodel$tpms, pop_mat = epimodel$pop_mat, eigen_vals = epimodel$eigen_values, eigen_vecs = epimodel$eigen_vectors, inverse_vecs = epimodel$inv_eigen_vectors, irm_keys = epimodel$keys)
                              
                              # TPM product subsequences
                              tpmProdSeqs(tpm_prods = epimodel$tpm_products, tpms = epimodel$tpms, obs_time_inds = epimodel$obs_time_inds)
                              
                              # Emission matrix
                              epimodel$emission_mat <- build_emission_mat(emission_mat = epimodel$emission_mat, epimodel = epimodel)
                              
                              # FB matrices
                              buildFBMats(epimodel$fb_mats, epimodel$tpm_products, epimodel$emission_mat, epimodel$initdist, epimodel$obs_time_inds)                        
                              # set variable for where to start storing additional state changes
                              epimodel$subj_row_ind <- epimodel$ind_final_config + 1
                              
                              # sample status at observation times
                              epimodel$config_mat[, subjects[j]] <- sample_at_obs_times(path = epimodel$config_mat[,subjects[j]], epimodel = epimodel)
                              
                              # sample status at event times
                              epimodel$config_mat[, subjects[j]] <- sample_at_event_times(path = epimodel$config_mat[,subjects[j]], epimodel = epimodel) 
                              
                              # sample paths in inter-event intervals
                              epimodel <- sample_path(epimodel = epimodel, subject = subjects[j])
                              
                              # insert the proposed trajectory
                              epimodel <- insert_trajectory(epimodel = epimodel, subject = subjects[j], reinsertion = FALSE)
                              
                              # accept or reject the proposed path
                              epimodel <- MH_accept_reject(epimodel = epimodel, subject = subjects[j], iter = j)
                              
                              epimodel$keys <- retrieveKeys(1:epimodel$ind_final_config, epimodel$irm_key_lookup, epimodel$pop_mat, epimodel$index_state_num)
                              
                              expect_equal(calc_pop_likelihood(epimodel, log = TRUE), epimodel$likelihoods$pop_likelihood_cur)
                    }
                    
          }
          
})