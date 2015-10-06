#' Update the observation matrix and the indices of observation times.
#' 
#' @inheritParams fit_epimodel
#'   
#' @return updated observation matrix and obs_time_inds vector in the epimodel
#'   environment
#' @export

update_obs_mat <- function(epimodel) {
          
          # set the vector of observation time indices
          epimodel$.obs_time_inds <- which(epimodel$config_mat[,"ID"] == 0)
          
          # update the compartment counts in the observation matrix
          epimodel$obs_mat[, epimodel$meas_vars_aug] <- epimodel$config_mat[epimodel$.obs_time_inds, epimodel$meas_vars]
          
}