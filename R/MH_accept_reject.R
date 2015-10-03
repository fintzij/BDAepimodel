#' Accept or reject a proposed subject-level trajectory via MH and update the
#' appropriate objects in the epimodel environment
#' 
#' @inheritParams draw_trajec
#'   
#' @return updated epimodel environment

MH_accept_reject <- function(epimodel, subject, subj_ID) {
          
          if((epimodel$likelihoods$pop_likelihood_new == -Inf ) | (epimodel$likelihoods$subj_likelihood_new == -Inf)) {
                    
                    acceptance_prob <- -Inf
                    
          } else {
                    
                    acceptance_prob <- epimodel$likelihoods$pop_likelihood_new - epimodel$likelihoods$pop_likelihood_cur + epimodel$likelihoods$subj_likelihood_cur - epimodel$likelihoods$subj_likelihood_new
                    
          }
          
          # if we accept the trajectory, need to update the observation matrix
          # and set population likelihood to new value
          if(acceptance_prob > 0 || acceptance_prob > log(runif(1))) {
                    
                    # update observation matrix to reflect new trajectory
                    epimodel$obs_mat[,paste0(epimodel$meas_vars, "_augmented")] <- epimodel$config_mat[epimodel$.obs_time_inds, epimodel$meas_vars]
                    
                    # set population level likelihood to the new value
                    epimodel$likelihoods$pop_likelihood_cur <- epimodel$likelihoods$pop_likelihood_new
                    
                    
          } else {
                    # remove the new trajectory
                    remove_trajectory(epimodel, subject = subject)
                    
                    reinsert_path(epimodel, subject = subject, subj_ID = subj_ID)
                    
                    # update the observation matrix
                    epimodel$obs_mat[,paste0(epimodel$meas_vars, "_augmented")] <- epimodel$config_mat[epimodel$.obs_time_inds, epimodel$meas_vars]
                    
          }
}