#' Initialize a set of trajectories if an initial set of paths not provided.
#'
#' @inheritParams draw_trajec
#'
#' @return epimodel object with initialized configuration matrix
#' @export
#'

init_trajectories <- function(epimodel) {
          
          if(is.null(epimodel$dat)) {
                    stop(sQuote("dat"), "must be specified within the epimodel")
          }
          
          epimod <- epimodel
          epimodel <- simulate_epimodel(epimodel = epimod, 
                                        obstimes = epimod$obstimes,
                                        meas_vars = epimod$meas_vars,
                                        r_meas_process = epimod$r_meas_process,
                                        return_config = TRUE,
                                        trim = FALSE,
                                        lump = TRUE)
          
          epimodel$dat <- epimod$dat
          epimodel$obs_mat[, paste0(epimodel$meas_vars, "_observed")] <- epimodel$dat[,epimodel$meas_vars]
          
          if(any(epimodel$d_meas_process(state = epimodel$obs_mat, meas_vars = epimodel$meas_vars, params = epimodel$params, log = FALSE) == 0)) {
                    while(any(epimodel$d_meas_process(state = epimodel$obs_mat, meas_vars = epimodel$meas_vars, params = epimodel$params, log = FALSE) == 0)) {
                              
                              epimod <- epimodel
                              epimodel <- simulate_epimodel(epimodel = epimod, 
                                                            obstimes = epimod$obstimes,
                                                            meas_vars = epimod$meas_vars,
                                                            r_meas_process = epimod$r_meas_process,
                                                            return_config = TRUE,
                                                            trim = FALSE,
                                                            lump = TRUE)
                              
                              epimodel$dat <- epimod$dat
                              epimodel$obs_mat[, paste0(epimodel$meas_vars, "_observed")] <- epimodel$dat[,epimodel$meas_vars]
                              
                    }
          }
          
          return(epimodel)
}