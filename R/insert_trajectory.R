#' Insert a newly sampled subject-level trajectory so its contribution is
#' reflected in the configuration matrix.
#' 
#' This function also updates several bookkeeping objects - .ind_final_config
#' 
#' @inheritParams sample_path
#'   
#' @return updated configuration matrix in the epimodel environment

insert_trajectory <- function(epimodel, subject, subj_ID) {
          
          # first re-order the configuration matrix if needed, then copy the
          # preceding rows for the new trajectory
          if(epimodel$.ind_final_config != (epimodel$.subj_row_ind - 1)) {
                    
                    # reorder the matrix
                    epimodel$config_mat[1:(epimodel$.subj_row_ind - 1),] <- epimodel$config_mat[order(epimodel$config_mat[1:(epimodel$.subj_row_ind - 1), "time"]), ]
                    
                    # copy the states from the previous event times
                    epimodel$config_mat[which(epimodel$config_mat[,"ID"] == subject), epimodel$states] <- epimodel$config_mat[which(epimodel$config_mat[,"ID"] == subject) - 1, epimodel$states]
                    
                    # copy the configurations from the previous event times
                    epimodel$config_mat[which(epimodel$config_mat[,"ID"] == subject), epimodel$.config_inds[-subject]] <- epimodel$config_mat[which(epimodel$config_mat[,"ID"] == subject) - 1, epimodel$.config_inds[-subject]] 
                    
                    # set the index for the final configuration
                    epimodel$.ind_final_config <- epimodel$.subj_row_ind - 1
                    
                    # set the vector of observation time indices
                    epimodel$.obs_time_inds <- which(epimodel$config_mat[,"ID"] == 0)
          }
          
          # add the subject's contribution back into the compartment counts
          .rows_to_update <- 1:epimodel$.ind_final_config
          
          epimodel$config_mat[, epimodel$states][cbind(.rows_to_update, epimodel$config_mat[, subj_ID][.rows_to_update])] <- epimodel$config_mat[,epimodel$states][cbind(.rows_to_update, epimodel$config_mat[, subj_ID][.rows_to_update])] + 1
}