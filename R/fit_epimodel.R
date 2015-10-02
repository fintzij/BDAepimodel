#' Main wrapper to fit a stochastic epimodel using the Bayesian data
#' augmentation algorithm.
#' 
#' @inheritParams simulate_epimodel
#'   
#' @return list containing the epimodel object and a results list. 
#' @export
#' 
fit_epimodel <- function(epimodel) {
          
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
          
          # generate .niter parameter samples
          for(k in 1:.niter) {
                    
                    # re-compute the array of rate matrices - only needs to be
                    # recomputed every time parameters are updated.
                    build_irm(.epimodel)
                    
                    # reset the vector of tpms to build (set all to false at the end of the build), then build the tpm
                    # sequences
                    .epimodel$.tpms_to_build <- get_tpms_to_build(.epimodel)
                    build_tpm_seqs(.epimodel)

                    # re-compute the population level likelihood, measurement
                    # process likelihood, and complete data log-likelihood
                    epimodel$likelihoods$pop_likelihood_cur <- calc_pop_likelihood(epimodel = .epimodel, log = TRUE)
                    epimodel$likelihoods$obs_likelihood_cur <- calc_obs_likelihood(epimodel = .epimodel, log = TRUE)

                    # choose which subjects should be redrawn
                    .subjects <- sample.int(n = .epimodel$popsize, size = .configs_to_redraw, replace = .config_replacement)

                    # cycle through subject-level trajectories to be re-drawn
                    for(j in 1:.configs_to_redraw) {
                              
                              # update_tpms to reflect removal of subject from
                              # the configuration
                              update_tpms(.epimodel, subject = .subjects[j], direction = "removal")
                              
                              # remove trajectory from the counts in config_mat 
                              # and obs_mat, and update the tpm sequences to 
                              # reflect the removal. 
                              remove_trajectory(.epimodel, subject = .subjects[j])
                              
                              # update instatiate missing IRMs, update the
                              # emission probability and FB matrices.
                              update_matrices(.epimodel, subject = .subjects[j])
                              
                              # draw a new subject level trajectory
                              draw_trajec(.epimodel, subject = .subjects[j])
                              
                              # check the irms to ensure that no additional
                              # matrices and decompositions are required
                              check_irm(.epimodel)
                              
                              # update the tpms to reflect the insertion of a
                              # subject-level trajectory into the configuration
                              update_tpms(.epimodel, subject = .subjects[j], direction = "insertion")
                              
                    }
                    
          }
          
          results <- init_results(epimodel)
}