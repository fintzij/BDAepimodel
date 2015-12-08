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
          
          epimod <- simulate_epimodel(epimodel = epimod, lump = TRUE, trim = FALSE)
          
          epimod$obs_mat[,epimodel$meas_vars_obs] <- epimod$dat[,epimodel$meas_vars] <- epimodel$dat[,epimodel$meas_vars]
          
          if(!all(epimod$d_meas_process(state = epimod$obs_mat, meas_vars = epimod$meas_vars, params = epimod$params, log = FALSE) > 0)) {
                    while(!all(epimod$d_meas_process(state = epimod$obs_mat, meas_vars = epimod$meas_vars, params = epimod$params, log = FALSE) > 0)) {
                              epimod <- epimodel
                              epimod <- simulate_epimodel(epimodel = epimod, lump = TRUE, trim = TRUE)
                              epimod$obs_mat[,epimodel$meas_vars_obs] <- epimod$dat[,epimod$meas_vars] <- epimodel$dat[,epimodel$meas_vars]
                    }
          }
          
          epimodel <- epimod
          
          return(epimodel)
}