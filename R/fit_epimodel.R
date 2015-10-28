#' Main wrapper to fit a stochastic epimodel using the Bayesian data
#' augmentation algorithm.
#'
#' @inheritParams simulate_epimodel
#' @param monitor TRUE/FALSE for whether to print the iteration number and
#'   produce a log-likelihood plot every time a trajectory is saved.
#'
#' @return list containing the epimodel object and a results list.
#' @export
#'
fit_epimodel <- function(epimodel, monitor = FALSE) {
        
        compiler::compilePKGS(enable=TRUE)
        compiler::enableJIT(3)
        
        # check that the simulation settings have been set
        if(is.null(epimodel$sim_settings)) {
                stop(sQuote("sim_settings"), "must be specified in the epimodel object.")
        }
        
        if(is.null(epimodel$meas_vars)) {
                stop(sQuote("meas_vars"), "must be specified within the epimodel object.")
        }
        
        if(is.null(epimodel$config_mat)) {
                warning("An initial configuration was not provided so one has been generated.")
                epimodel <- init_augmentation(epimodel)
        }
        
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
                
                # compute the rates
                rates <- do.call(cbind, lapply(epimodel$rates, do.call, list(state = epimodel$irm_key_lookup, params = epimodel$params)))
                
                # update the rate matrices
                buildRateArray(irm_array = epimodel$irm, rates = rates, flow_inds = epimodel$flow_inds)
                
                # update the eigen decompositions
                buildEigenArray(eigenvals = epimodel$eigen_values, eigenvecs = epimodel$eigen_vectors, inversevecs = epimodel$inv_eigen_vectors, irm_array = epimodel$irm)
                
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
                        epimodel$fb_mats <- build_fb_mats(fb_mats = epimodel$fb_mats, epimodel = epimodel)
                        
                        # draw a new subject level trajectory
                        draw_trajec(epimodel, subject = subjects[j], iter = j)
                        
                }
                
                # record the acceptance rate for trajectories
                results$accepts[k - 1,"paths"] <- mean(epimodel$path_accept_vec)
                
                # update the measurement process likelihood
                epimodel$likelihoods$obs_likelihood <- calc_obs_likelihood(epimodel, log = TRUE)
                
                # get the current values of the parameters
                .params_cur <- epimodel$params
                
                # sample new parameters
                for(t in 1:length(kernel)) {
                        kernel[[t]](epimodel)
                }
                
                # sample new values for the initial distribution parameters
                epimodel$initdist_kernel(epimodel)
                
                # record acceptances/rejections for the parameter updates
                results$accepts[k - 1, names(epimodel$params)] <- as.numeric(epimodel$params != .params_cur)
                
                # save parameter values and log-likelihood
                if(k %% save_params_every == 0) {
                        
                        results$params[k %/% save_params_every, ] <- epimodel$params
                        results$log_likelihood[k %/% save_params_every] <- sum(epimodel$likelihoods$obs_likelihood, epimodel$likelihoods$pop_likelihood_cur)
                        
                }
                
                # save the configuration matrix
                if(k %% save_configs_every == 0) {
                        
                        results$configs[[k %/% save_configs_every]] <- epimodel$config_mat[complete.cases(epimodel$config_mat),]
                }
                
                if(monitor && (k%%save_configs_every) == 0) {
                        print(k)
                        # ts.plot(results$log_likelihood, xlab = "Iteration", ylab = "log-likelihood")
                }
        }
        
        end.time <- Sys.time()
        
        results$time <- difftime(end.time, start.time, units = "hours")
        
        epimodel <- list(model = list(dat = epimodel$dat,
                                      time_var = epimodel$time_var,
                                      meas_vars = epimodel$meas_vars,
                                      obstimes = epimodel$obstimes,
                                      popsize = epimodel$popsize,
                                      states = epimodel$states,
                                      params = apply(results$params, 2, median)),
                         settings = epimodel$sim_settings,
                         results = results)
        
        # clean up internal objects
        rm(epimodel)
        gc()
        
        return(epimodel)
}
