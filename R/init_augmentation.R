#' Generate an initial configuration that is concordant with the data, updating
#' the configuration and observation matrices.
#' 
#' @param epimodel
#'   
#' @return initialized epimodel object
#' @export
#' 
init_augmentation <- function(epimodel) {
          
          epimod <- epimodel
          
          epimodel <- simulate_epimodel(epimodel = epimod, lump = TRUE, trim = FALSE)
          
          epimodel$obs_mat[,epimodel$meas_vars_obs] <- epimodel$dat[,epimodel$meas_vars]
          
          if(!all(epimodel$d_meas_process(state = epimodel$obs_mat, meas_vars = epimodel$meas_vars, params = epimodel$params, log = FALSE) > 0)) {
                    while(!all(epimodel$d_meas_process(state = epimodel$obs_mat, meas_vars = epimodel$meas_vars, params = epimodel$params, log = FALSE) > 0)) {
                              epimodel <- epimod
                              epimodel <- simulate_epimodel(epimodel = epimod, lump = TRUE, trim = TRUE)
                    }
          }
          
          return(epimodel)
}