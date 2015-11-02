#' Accept or reject a proposed subject-level trajectory via MH and update the
#' appropriate objects in the epimodel environment
#' 
#' @inheritParams draw_trajec
#'   
#' @return updated epimodel environment
#' @export

MH_accept_reject <- function(epimodel, subject, iter) {
          
          if((epimodel$likelihoods$pop_likelihood_new == -Inf ) & (epimodel$likelihoods$subj_likelihood_new == -Inf)) {
                    
                    acceptance_prob <- -Inf
                    
          } else {
                    
                    acceptance_prob <- epimodel$likelihoods$pop_likelihood_new - epimodel$likelihoods$pop_likelihood_cur + epimodel$likelihoods$subj_likelihood_cur - epimodel$likelihoods$subj_likelihood_new
                    
          }
          
          # if we accept the trajectory, need to update the observation matrix
          # and set population likelihood to new value
          if(acceptance_prob > 0 || acceptance_prob > log(runif(1))) {
                    
                    # set population level likelihood to the new value
                    epimodel$likelihoods$pop_likelihood_cur <- epimodel$likelihoods$pop_likelihood_new
                    
                    # record acceptance
                    epimodel$path_accept_vec[iter] <- 1
                    
          } else {
                    # the rejected path is being re-inserted into the config mat since it is stored in remove trajectory. need to save it separately first
                    # remove the proposed path, not overwriting the current path
                    epimodel <- remove_trajectory(epimodel, subject = subject, save_path = FALSE)
                    
                    # reinsert the previous path
                    epimodel <- reinsert_path(epimodel, subject = subject)
                    
                    # record rejection
                    epimodel$path_accept_vec[iter] <- 0
          }
          
          # update counts in the obervation matrix
          epimodel$obs_mat[, epimodel$meas_vars_aug] <- epimodel$pop_mat[epimodel$obs_time_inds, epimodel$meas_vars]
          
          # update the initial configuration vector
          epimodel$init_config[subject] = epimodel$subj_path[1]
          
          return(epimodel)
}