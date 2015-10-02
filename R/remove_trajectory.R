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
          .subj_ID         <- paste0(".X", subject)
          .subj_inds          <- which(epimodel$config_mat[,"ID"] == subject)
          
          epimodel$.path_cur  <- epimodel$config_mat[c(1, .subj_inds, epimodel$.ind_final_config), c("time", "ID", "Event", .subj_ID)]
          
          # remove the contribution to the compartment counts in the configuration matrix
          .rows_to_update <- 1:epimodel$.ind_final_config
          
          epimodel$config_mat[, epimodel$states][cbind(.rows_to_update, epimodel$config_mat[, .subj_ID][.rows_to_update])] <- epimodel$config_mat[,epimodel$states][cbind(.rows_to_update, epimodel$config_mat[, .subj_ID][.rows_to_update])] - 1
          
          # compute the likelihood of the subject's trajectory from a
          # time-inhomogeneous CTMC
          epimodel$likelihoods$subj_likelihood_cur <- calc_subj_likelihood(epimodel = epimodel, subject = subject, subj_ID = .subj_ID, log = TRUE)
          
          # remove the rows relating to the trajectory from the configuration matrix
          epimodel$config_mat[.subj_inds,] <- NA
          epimodel$config_mat <- epimodel$config_mat[order(epimodel$config_mat[,"time"]),]
          
          # set the index of the final configuration
          epimodel$.ind_final_config <- epimodel$.ind_final_config - (nrow(epimodel$.path_cur) - 2)
          
          # set the vector of observation time indices
          epimodel$.obs_time_inds <- which(epimodel$config_mat[,"ID"] == 0)
          
          # remove the contribution to the compartment counts in the observation matrix
          epimodel$obs_mat[, paste0(epimodel$meas_vars, "_augmented")] <- (epimodel$config_mat[, epimodel$meas_vars][1:epimodel$.ind_final_config])[epimodel$.obs_time_inds]
          
}