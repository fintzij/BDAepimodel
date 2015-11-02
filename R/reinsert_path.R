#' Re-insert the current path back into the configuration matrix.
#' 
#' Note that this gets called after the proposed path has been removed.
#' Therefore, .ind_final_config is already updated.
#' 
#' @inheritParams draw_trajec
#'   
#' @return updated configuration matrix
#' 
#' @export

reinsert_path <- function(epimodel, subject) {
          
          # find intervals in .path_cur times for matching indices
          inds <- findInterval(epimodel$pop_mat[1:epimodel$ind_final_config,"time"], epimodel$path_cur[,1])
          
          # copy the state from the current path into the config_mat
          epimodel$subj_path[1:epimodel$ind_final_config] <- as.integer(epimodel$path_cur[inds, 3])
          
          # set the number of jumps to the number in the original trajectory
          epimodel$n_jumps <- nrow(epimodel$path_cur) - 2
          
          # if there were transitions in the path, copy them to the
          # configuration matrix and call insert_trajectory
          if(epimodel$n_jumps > 0) {
                    
                    insertPath(epimodel$path_cur[2:(nrow(epimodel$path_cur) - 1), , drop = FALSE], subject, epimodel$pop_mat, epimodel$subj_path, epimodel$ind_final_config)

          }
          
          # insert_trajectory will insert the transitions and reorder the matrix
          epimodel <- insert_trajectory(epimodel, subject = subject, reinsertion = TRUE)
          
          return(epimodel)
}