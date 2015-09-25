#' Remove the contribution of a subject from the compartment counts in the 
#' \code{config_mat} and \code{obs_mat} objects.
#' 
#' @inheritParams simulate_epimodel
#' @param subject ID for subject to be removed
#'   
#' @return \code{config_mat} and \code{obs_mat} objects within the .epimodel
#'   environment.

remove_trajectory <- function(epimodel, subject) {
          
          .rows_to_replace <- 1:epimodel$.ind_final_config
  
          # remove the contribution to the compartment counts in the configuration matrix
          epimodel$config_mat[, epimodel$states][cbind(.rows_to_replace, epimodel$config_mat[, subject][.rows_to_replace])] <- epimodel$config_mat[, epimodel$states][cbind(.rows_to_replace, epimodel$config_mat[, subject][.rows_to_replace])] - 1
          
          # remove the contribution to the compartment counts in the observation matrix
          epimodel$obs_mat[, paste0(epimodel$meas_vars, "_augmented")] <- (epimodel$config_mat[, epimodel$meas_vars][.rows_to_replace])[epimodel$config_mat[,"Event"][.rows_to_replace] == 0]

}