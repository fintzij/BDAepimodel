#' Re-insert the current path back into the configuration matrix.
#' 
#' Note that this gets called after the proposed path has been removed.
#' Therefore, .ind_final_config is already updated.
#' 
#' @inheritParams draw_trajec
#'   
#' @return updated configuration matrix

reinsert_path <- function(epimodel, subject, subj_ID) {
          
          # find intervals in .path_cur times for matching indices
          .inds <- findInterval(epimodel$config_mat[1:epimodel$.ind_final_config,"time"], epimodel$.path_cur[,"time"])
          
          # copy the state from the current path into the config_mat
          epimodel$config_mat[1:epimodel$.ind_final_config, subj_ID] <- epimodel$.path_cur[.inds, subj_ID]
          
          # if there were transitions in the path, copy them to the
          # configuration matrix and call insert_trajectory
          if(nrow(epimodel$.path_cur) > 2) {
                    
                    epimodel$config_mat[(epimodel$.ind_final_config + 1) : ((epimodel$.ind_final_config + nrow(epimodel$.path_cur) - 2)), c("time", "ID", "Event", subj_ID)] <- epimodel$.path_cur[2:(nrow(epimodel$.path_cur) - 1), ]
                    
                    # set the final index of the subject path in the config_mat
                    epimodel$.subj_row_ind <- epimodel$.ind_final_config + nrow(epimodel$.path_cur) - 1
                    
                    # insert_trajectory will insert the transitions and reorder the matrix
                    insert_trajectory(epimodel, subject = subject, subj_ID = subj_ID)
                    
          }
}