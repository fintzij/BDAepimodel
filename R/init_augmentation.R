#' Generate an initial configuration that is concordant with the data, updating 
#' the configuration and observation matrices.
#' 
#' @param epimodel
#'   
#' @return initialized epimodel object
#' @export
#' 
init_augmentation <- function(epimodel) {
          
          SIR  <- identical(epimodel$states, c("S","I","R")) && (nrow(epimodel$flow) == 2)

          if(is.null(epimodel$sim_settings$init_popsize)) {
                    init_popsize <- epimodel$popsize
          } else {
                    init_popsize <- epimodel$sim_settings$init_popsize
          }
          
          epimod <- epimodel
          
          if(SIR) {
                    init_config <- simulate_SIR(epimodel$obstimes, epimodel$params, init_popsize, trim = FALSE)        
                    epimod$obs_mat <- init_config$obs_mat
                    epimod$pop_mat <- init_config$pop_mat
                    epimod$init_config <- rep(1:epimodel$num_states, epimod$pop_mat[1,epimod$states])
          } else {
                    epimod <- suppressWarnings(simulate_epimodel(epimodel = epimod, lump = TRUE, trim = FALSE))
          }

          epimod$obs_mat[,epimodel$meas_vars_obs] <- epimod$dat[,epimodel$meas_vars] <- epimodel$dat[,epimodel$meas_vars]
          
          initialization_attempts <- 1
          
          # ensure all likelihood values are positive and there are no nan values
          liks <- suppressWarnings(epimod$d_meas_process(state = epimod$obs_mat, 
                                                         meas_vars = epimod$meas_vars, 
                                                         params = epimod$params, log = FALSE))
          valid_init <- all(liks > 0) && !any(is.nan(liks))
          
          if(!valid_init) {
                    while(!valid_init) {
                              
                              initialization_attempts <- initialization_attempts + 1
                              
                              if(initialization_attempts%%10 == 0) {
                              print(paste0("Initialization attempt ", initialization_attempts))
                              }
                              
                              if(initialization_attempts > 1000) {
                                        stop("Initialization failed. Try other parameter values.")
                              }
                              
                              epimod <- epimodel
                              
                              if(SIR) {
                                        init_config <- simulate_SIR(epimodel$obstimes, epimodel$params, init_popsize, trim = FALSE)        
                                        epimod$obs_mat <- init_config$obs_mat
                                        epimod$pop_mat <- init_config$pop_mat
                                        epimod$init_config <- rep(1:epimodel$num_states, epimod$pop_mat[1,epimod$states])
                              } else {
                                        epimod <- suppressWarnings(simulate_epimodel(epimodel = epimod, lump = TRUE, trim = FALSE))
                              }
                              
                              epimod$obs_mat[,epimodel$meas_vars_obs] <- epimod$dat[,epimod$meas_vars] <- epimodel$dat[,epimodel$meas_vars]
                              
                              liks <- suppressWarnings(epimod$d_meas_process(state = epimod$obs_mat, 
                                                                             meas_vars = epimod$meas_vars, 
                                                                             params = epimod$params, log = FALSE))
                              valid_init <- all(liks > 0) && !any(is.nan(liks))
                    }
          }
          
          epimodel <- epimod
          
          if(init_popsize < epimodel$popsize) {
                    state_probs <- epimodel$sim_settings$compartment_dist[epimodel$states]
                    other_subj_inits <-
                              replicate(
                                        n = epimodel$popsize - epimodel$sim_settings$init_popsize,
                                        expr = sample.int(epimodel$num_states, 1, prob = state_probs)
                              )
                    other_subj_statecounts <-
                              sapply(1:epimodel$num_states, function(x)
                                        sum(other_subj_inits == x))
                    epimodel$init_config <- c(epimodel$init_config, other_subj_inits)
                    epimodel$pop_mat[, 4:ncol(epimodel$pop_mat)] <-
                              sweep(epimodel$pop_mat[, 4:ncol(epimodel$pop_mat)], 2, other_subj_statecounts, FUN =
                                              "+")
          }
          
          return(epimodel)
}