#' Remove the contribution of a subject from the compartment counts in the 
#' \code{config_mat} and \code{obs_mat} objects.
#' 
#' @inheritParams simulate_epimodel
#' @param subject ID for subject to be removed
#'   
#' @return \code{config_mat} and \code{obs_mat} objects within the .epimodel
#'   environment.

remove_trajectory <- function(epimodel, subject) {
          
          # extract the current trajectory
          .subject_ID         <- paste0(".X", subject)
          .subj_inds          <- which(epimodel$config_mat[,"ID"] == subject)
          
          epimodel$.path_cur  <- epimodel$config_mat[c(1, .subj_inds, epimodel$.ind_final_config), c("time", "ID", "Event", epimodel$states, .subject_ID)]
          
          # remove the rows relating to the trajectory from the configuration matrix
          epimodel$config_mat[.subj_inds,] <- NA
          epimodel$config_mat <- epimodel$config_mat[order(epimodel$config_mat[,"time"]),]
          
          # set the index of the final configuration
          epimodel$.ind_final_config <- epimodel$.ind_final_config - (nrow(epimodel$.path_cur) - 2)
          .rows_to_update <- 1:epimodel$.ind_final_config
          
          # set the vector of observation time indices
          epimodel$.obs_time_inds <- which(epimodel$config_mat[,"ID"] == 0)
  
          # remove the contribution to the compartment counts in the configuration matrix
          epimodel$config_mat[, epimodel$states][cbind(.rows_to_update, epimodel$config_mat[, .subject_ID][.rows_to_update])] <- epimodel$config_mat[,epimodel$states][cbind(.rows_to_update, epimodel$config_mat[, .subject_ID][.rows_to_update])] - 1
          
          # remove the contribution to the compartment counts in the observation matrix
          epimodel$obs_mat[, paste0(epimodel$meas_vars, "_augmented")] <- (epimodel$config_mat[, epimodel$meas_vars][.rows_to_update])[epimodel$.obs_time_inds]
          
          # re-order the .tpms and .tpm_products objects, setting the offending
          # matrices to null and placing them at the end
          .tpm_order <- c(setdiff(1:length(epimodel$.tpms), .subj_inds), .subj_inds)
          epimodel$.tpms[.subj_inds] <- epimodel$.tpm_products[.subj_inds] <- list(NULL)
          epimodel$.tpms <- epimodel$.tpms[.tpm_order]
          epimodel$.tpm_products <- epimodel$.tpm_products[.tpm_order]
}