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
          
          suppressMessages(require(MASS, quietly = TRUE))
          suppressMessages(require(MCMCpack, quietly = TRUE))
          suppressMessages(require(compiler, quietly = TRUE))
          
          compilePKGS(enable=TRUE)
          enableJIT(3)
          
          # check that the simulation settings have been set
          if(is.null(epimodel$sim_settings)) {
                    stop(sQuote("sim_settings"), "must be specified in the epimodel object.")
          }
          
          if(is.null(epimodel$meas_vars)) {
                    stop(sQuote("meas_vars"), "must be specified within the epimodel object.")
          }
          
          if(is.null(epimodel$config_mat)) {
                    warning("An initial configuration was not provided so one has been generated.")
                    #############################################
                    ### NEED TO WRITE INITIALIZATION FUNCTION ###
                    #############################################
                    epimodel$config_mat
          }
          
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
          
          start.time = Sys.time()
          # generate .niter parameter samples
          for(k in 2:niter) {
                    if(k%%5 == 0) print(k)
                    # re-compute the array of rate matrices - only needs to be
                    # recomputed every time parameters are updated.
                    build_irm(.epimodel)

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
                    
                    # update the measurement process likelihood
                    .epimodel$likelihoods$obs_likelihood <- calc_obs_likelihood(.epimodel, log = TRUE)
                    
                    # sample new parameters
                    for(t in 1:length(kernel)) {
                              kernel[[t]](.epimodel)
                    }
                    
                    # sample new values for the initial distribution parameters
                    .epimodel$initdist_kernel(.epimodel)
                    
                    # save parameter values and log-likelihood
                    if(k %% save_params_every == 0) {
                              
                              results$params[k %/% save_params_every, ] <- .epimodel$params
                              results$log_likelihood[k %/% save_params_every] <- sum(.epimodel$likelihoods$obs_likelihood, .epimodel$likelihoods$pop_likelihood_cur)
                              
                    }
                    
                    # save the configuration matrix
                    if(k %% save_configs_every == 0) {
                              
                              results$configs[[k %/% save_configs_every]] <- .epimodel$config_mat[complete.cases(.epimodel$config_mat),]
                    }
                    
                    if(monitor && (k%%save_configs_every) == 0) {
                              print(k)
                              ts.plot(results$log_likelihood, xlab = "Iteration", ylab = "log-likelihood")
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
                           settings = .epimodel$sim_settings,
                           results = results)
          
          # clean up internal objects
          rm(.epimodel)
          gc()
          enableJIT(0)
          
          return(epimodel)
}