#' Generate an initial configuration that is concordant with the data, updating 
#' the configuration and observation matrices.
#' 
#' @param epimodel
#'   
#' @return initialized epimodel object
#' @export
#' 
init_augmentation <- function(epimodel) {
          
          SIR <- all(epimodel$states == c("S","I","R"))
          
          epimod <- epimodel
          
          if(SIR) {
                    init_config <- simulate_SIR(epimodel$obstimes, epimodel$params, epimodel$popsize, trim = FALSE)        
                    epimod$obs_mat <- init_config$obs_mat
                    epimod$pop_mat <- init_config$pop_mat
                    epimod$init_config <- rep(1:epimodel$num_states, epimod$pop_mat[1,epimod$states])
          } else {
                    epimod <- simulate_epimodel(epimodel = epimod, lump = TRUE, trim = FALSE)
          }         

          epimod$obs_mat[,epimodel$meas_vars_obs] <- epimod$dat[,epimodel$meas_vars] <- epimodel$dat[,epimodel$meas_vars]
          
          if(!all(epimod$d_meas_process(state = epimod$obs_mat, meas_vars = epimod$meas_vars, params = epimod$params, log = FALSE) > 0)) {
                    while(!all(epimod$d_meas_process(state = epimod$obs_mat, meas_vars = epimod$meas_vars, params = epimod$params, log = FALSE) > 0)) {
                              
                              epimod <- epimodel
                              
                              if(SIR) {
                                        init_config <- simulate_SIR(epimodel$obstimes, epimodel$params, epimodel$popsize, trim = FALSE)        
                                        epimod$obs_mat <- init_config$obs_mat
                                        epimod$pop_mat <- init_config$pop_mat
                                        epimod$init_config <- rep(1:epimodel$num_states, epimod$pop_mat[1,epimod$states])
                              } else {
                                        epimod <- simulate_epimodel(epimodel = epimod, lump = TRUE, trim = FALSE)
                              }
                              epimod$obs_mat[,epimodel$meas_vars_obs] <- epimod$dat[,epimod$meas_vars] <- epimodel$dat[,epimodel$meas_vars]
                    }
          }
          
          epimodel <- epimod
          
          return(epimodel)
}