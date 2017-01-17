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

          # check that the simulation settings have been set
          if(is.null(epimodel$sim_settings)) {
                stop(sQuote("sim_settings"), "must be specified in the epimodel object.")
          }

          if(is.null(epimodel$meas_vars)) {
                stop(sQuote("meas_vars"), "must be specified within the epimodel object.")
          }

          if(is.null(epimodel$pop_mat)) {
                    epimodel <- init_augmentation(epimodel)
                    print("Configuration initialized. Beginning MCMC.")
          } else {
                    # generate the initial configuration of subject labels at time t0
                    epimodel$init_config <- rep(1:epimodel$num_states, epimodel$pop_mat[1,epimodel$states])

                    # initialize the observation matrix
                    epimodel$obs_mat <- init_obs_mat(epimodel)
          }

          if(!is.null(epimodel$sim_settings$post_init_params)) {
                    epimodel$params <- epimodel$sim_settings$post_init_params
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
          results               <- init_results(epimodel)

          # set seed
          results$seed          <- epimodel$sim_settings$seed
          set.seed(results$seed)

          # set the maximum run time
          if(is.null(epimodel$sim_settings$time_limit)) {
                  max_runtime <- as.difftime(Inf, units = "hours")
          } else {
                  max_runtime <- epimodel$sim_settings$time_limit
          }

          # store the initial values in the results matrix
          results$params[1, ]       <- epimodel$params
          results$log_likelihood[1] <- epimodel$likelihoods$obs_likelihood + epimodel$likelihoods$pop_likelihood_cur
          results$configs[[1]]      <- epimodel$pop_mat[1:epimodel$ind_final_config,]

          # initialize the vector for path acceptances
          epimodel$path_accept_vec <- rep(0, configs_to_redraw)

          # build the vector of probabilities for selecting which subjects to sample
          subj_probs <- subject_probabilities(epimodel)
          subjects   <- integer(configs_to_redraw)

          # get start time
          start.time = Sys.time()

          # generate niter parameter samples
          for(k in 2:niter) {

                    # choose which subjects should be redrawn
                    if(subj_probs$preferential_sampling) {

                              pref_group    <- as.logical(rbinom(configs_to_redraw, 1, subj_probs$pref_prob))
                              n_pref        <- sum(pref_group) # number of subjects from the preferential group

                              # sample IDs of subjects in the preferential and non-preferential groups
                              subjects_pref    <- sample(x = subj_probs$pref_IDs,
                                                         size = n_pref,
                                                         replace = config_replacement)
                              subjects_nonpref <- sample(x = subj_probs$nonpref_IDs,
                                                         size = configs_to_redraw - n_pref,
                                                         replace = config_replacement)

                              # fill out the vector of subject IDs
                              subjects[pref_group]  <- subjects_pref
                              subjects[!pref_group] <- subjects_nonpref

                    } else {
                              subjects <- sample.int(
                                        n       = epimodel$popsize,
                                        size    = configs_to_redraw,
                                        replace = config_replacement,
                                        prob    = subj_probs$subj_probs
                              )
                    }

                    # cycle through subject-level trajectories to be re-drawn
                    for(j in 1:configs_to_redraw) {

                              # remove trajectory from the counts in pop_mat
                              # and obs_mat, and update the tpm sequences to
                              # reflect the removal.
                              epimodel <- remove_trajectory(epimodel, subject = subjects[j], save_path = TRUE)

                              # TPM sequence
                              tpmSeqs(
                                        tpms            = epimodel$tpms,
                                        pop_mat         = epimodel$pop_mat,
                                        real_eigen_vals = epimodel$real_eigen_values,
                                        imag_eigen_vals = epimodel$imag_eigen_values,
                                        eigen_vecs      = epimodel$eigen_vectors,
                                        inverse_vecs    = epimodel$inv_eigen_vectors,
                                        irm_keys        = epimodel$keys,
                                        n_real_eigs     = epimodel$n_real_eigs,
                                        irms            = epimodel$irm
                              )

                              # TPM product subsequences
                              tpmProdSeqs(
                                        tpm_prods     = epimodel$tpm_products,
                                        tpms          = epimodel$tpms,
                                        obs_time_inds = epimodel$obs_time_inds
                              )

                              # Emission matrix
                              epimodel$emission_mat <- build_emission_mat(emission_mat = epimodel$emission_mat, epimodel = epimodel)

                              # FB matrices
                              buildFBMats(
                                        epimodel$fb_mats,
                                        epimodel$tpm_products,
                                        epimodel$emission_mat,
                                        epimodel$initdist,
                                        epimodel$obs_time_inds
                              )

                              # sample status at observation times
                              epimodel$subj_path <- sample_at_obs_times(path = epimodel$subj_path, epimodel = epimodel)

                              # sample status at event times
                              epimodel$subj_path <- sample_at_event_times(path = epimodel$subj_path, epimodel = epimodel)

                              # sample paths in inter-event intervals
                              path               <- sample_path(epimodel)
                              epimodel$n_jumps   <- nrow(path)

                              # insert the path
                              if(epimodel$n_jumps) {

                                        # add rows to the population level bookkeeping matrix if necessary
                                        if(epimodel$ind_final_config + epimodel$n_jumps > nrow(epimodel$pop_mat)) {
                                                  epimodel <- expand_pop_mat(epimodel, buffer_size = epimodel$n_jumps)
                                        }

                                        # insert the path
                                        insertPath(path,
                                                   subjects[j],
                                                   epimodel$pop_mat,
                                                   epimodel$subj_path,
                                                   epimodel$ind_final_config)
                              }

                              # insert the proposed trajectory
                              epimodel <- insert_trajectory(epimodel    = epimodel,
                                                            subject     = subjects[j],
                                                            reinsertion = FALSE)

                              # accept or reject the proposed path
                              epimodel <- MH_accept_reject(epimodel = epimodel, subject = subjects[j], iter = j)
                    }

                    # record the acceptance rate for trajectories
                    results$accepts[k - 1,"paths"] <- mean(epimodel$path_accept_vec)

                    # update the measurement process likelihood
                    epimodel$likelihoods$obs_likelihood <- calc_obs_likelihood(epimodel, log = TRUE)

                    # get the current values of the parameters
                    epimodel$.params_cur <- epimodel$params
                    epimodel$.irm_cur    <- epimodel$irm

                    # retrieve the irm keys (last trajectory update could have
                    # been rejected, so keys would be incorrect)
                    epimodel$keys <-
                              retrieveKeys(
                                        1:epimodel$ind_final_config,
                                        epimodel$irm_key_lookup,
                                        epimodel$pop_mat,
                                        epimodel$index_state_num
                              )

                    # sample new parameters
                    for(t in 1:length(kernel)) {
                              epimodel <- kernel[[t]](epimodel)
                    }

                    # sample new values for the initial distribution parameters
                    epimodel$params     <- epimodel$initdist_kernel(epimodel)
                    epimodel$initdist   <- epimodel$params[epimodel$initdist_params]

                    # record acceptances/rejections for the parameter updates
                    results$accepts[k - 1, names(epimodel$params)] <- as.numeric(epimodel$params != epimodel$.params_cur)

                    # save parameter values and log-likelihood
                    if(k %% save_params_every == 0) {
                        results$params[k %/% save_params_every, ]       <- epimodel$params
                        results$log_likelihood[k %/% save_params_every] <- sum(epimodel$likelihoods$obs_likelihood,
                                                                               epimodel$likelihoods$pop_likelihood_cur)

                    }

                    # save the configuration matrix
                    if(k %% save_configs_every == 0) {
                        results$configs[[k %/% save_configs_every]] <- epimodel$pop_mat[complete.cases(epimodel$pop_mat),]
                    }

                    if(monitor && (k%%save_configs_every) == 0) {
                        print(k)
                        # ts.plot(results$log_likelihood, xlab = "Iteration", ylab = "log-likelihood")
                    }

                    if((k %% save_params_every == 0) && (difftime(Sys.time(), start.time) > max_runtime)) break
          }

          end.time <- Sys.time()

          results$time <- difftime(end.time, start.time, units = "hours")
          results$iterations_completed <- k

          epimodel <- list(dat = epimodel$dat,
                         time_var = epimodel$time_var,
                         meas_vars = epimodel$meas_vars,
                         obstimes = epimodel$obstimes,
                         popsize = epimodel$popsize,
                         states = epimodel$states,
                         params = apply(results$params, 2, median),
                         settings = epimodel$sim_settings,
                         results = results)

        return(epimodel)
}
