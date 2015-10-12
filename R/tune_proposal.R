#' Tune proposal distribution for randow walk M-H parameter proposals.
#' 
#' @param epimodel
#' @param epimodel_envir logical indicating whether the tuning is run within the
#'   epimodel environment
#'   
#' @return updated covariance matrix and scale factor within the epimodel
#'   environment, or a list with tuning results
#' @export
#' 

tune_proposal <- function(epimodel, epimodel_envir) {
          
          # set the number of iterations to be equal to the number of iterations
          # per tuning loop
          epimodel$sim_settings$niter <- epimodel$tune_settings$ntu
          
          # don't save trajectories
          epimodel$sim_settings$save_configs_every <- epimodel$tune_settings$ntu + 1
          
          # save parameters every iteration
          epimodel$sim_settings$save_params_every <- 1
          
          # initialize acceptance and parameter matrices
          accepts <- params <- matrix(0, nrow = epimodel$tune_settings$ntu, ncol = length(epimodel$tune_settings$tune_params))
          colnames(accepts) <- colnames(params) <- epimodel$tune_settings$tune_params
          
          # if not tuning within an epimodel environment, prepare an epimodel environment
          if(!is.environment(epimodel)) {
                    .epimodel <- prepare_epimodel(epimodel)
          }
          
          # run the tuning loops
          for(p in 1:.epimodel$tune_settings$maxtune) {
                    
                    print(paste("Tuning loop", p))
                    
                    # generate ntu parameter samples
                    for(k in 1:.epimodel$tune_settings$ntu) {
                              
                              # re-compute the array of rate matrices - only needs to be
                              # recomputed every time parameters are updated.
                              build_irm(.epimodel)
                              
                              # choose which subjects should be redrawn
                              .subjects <- sample.int(n = .epimodel$popsize, size = .epimodel$sim_settings$configs_to_redraw, replace = .epimodel$sim_settings$config_replacement)
                              
                              # cycle through subject-level trajectories to be re-drawn
                              for(j in 1:.epimodel$sim_settings$configs_to_redraw) {
                                        
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
                                        draw_trajec(.epimodel, subject = .subjects[j], iter = j) 
                                        
                              }
                              
                              # update the measurement process likelihood
                              .epimodel$likelihoods$obs_likelihood <- calc_obs_likelihood(.epimodel, log = TRUE)
                              
                              # get the current values of the parameters
                              .params_cur <- .epimodel$params
                              
                              # sample new parameters
                              for(t in 1:length(kernel)) {
                                        .epimodel$sim_settings$kernel[[t]](.epimodel)
                              }
                              
                              # sample new values for the initial distribution parameters
                              .epimodel$initdist_kernel(.epimodel)
                              
                              # record acceptances/rejections for the parameter updates
                              accepts[k, ] <- as.numeric(.epimodel$params[.epimodel$tune_settings$tune_params] == .params_cur[.epimodel$tune_settings$tune_params])
                              
                              # save the parameter values
                              if(is.null(.epimodel$sim_settings$to_estimation_scale)) {
                                        params[k, ] <- .epimodel$params[.epimodel$tune_settings$tune_params]
                                        
                              } else {
                                        params[k, ] <- .epimodel$sim_settings$to_estimation_scale(.epimodel$params, .epimodel)[.epimodel$tune_settings$tune_params]
                              }
                    }
                    
                    if(findInterval(mean(accepts), .epimodel$tune_settings$target_accept) & p >= .epimodel$tune_settings$mintune) break
                    
                    # compute the covariance matrix for the current run
                    if(.epimodel$tune_settings$tunewt != 0) {
                              
                              wt_cur <- .epimodel$tune_settings$tunewt / (1 + (p-1) * .epimodel$tune_settings$tunewt)
                              
                              .epimodel$sim_settings$cov_mtx <- wt_cur * cov(params) + (1 - wt_cur) * .epimodel$sim_settings$cov_mtx + diag(rnorm(length(.epimodel$tune_settings$tune_params), 0, 1e-6))
                    }
                    
                    # compute the scale factor
                    .epimodel$tune_settings$cov_scale <- .epimodel$tune_settings$cov_scale * qnorm(0.234 / 2) / qnorm(mean(accepts) / 2)
                    
          }
          
          
          if(epimodel_envir) {
                    .epimodel$sim_settings <- .epimodel$tune_settings$cov_scale^2 * .epimodel$tune_settings$cov_scale
          } else {
                    cov_scale <- .epimodel$tune_settings$cov_scale
                    cov_mtx <- .epimodel$sim_settings$cov_mtx
                    
                    rm(.epimodel)
                    gc()
                    
                    return(list(cov_mtx = cov_mtx,
                                cov_scale = cov_scale,
                                accept_rate = mean(accepts),
                                num_adaptations = p))
          }
          
} 