#' Remove the contribution of a subject from the compartment counts in the 
#' \code{config_mat} and \code{obs_mat} objects.
#' 
#' @inheritParams simulate_epimodel
#' @param subject ID for subject to be removed
#' @param save_path TRUE/FALSE for whether to save the current path in epimodel
#'   
#' @return \code{config_mat} and \code{obs_mat} objects within the .epimodel
#'   environment.
#'   
#' @export

remove_trajectory <- function(epimodel, subject, save_path) {
          
          # extract the current trajectory
          .subj_ID            <- paste0(".X", subject)
          .subj_inds          <- which(epimodel$config_mat[,"ID"] == subject)
          
          # save the current path if not removing it b/c of a M-H rejection
          if(save_path) {
                    epimodel$.path_cur  <- epimodel$config_mat[c(1, .subj_inds, epimodel$.ind_final_config), c("time", "ID", "Event", .subj_ID)]
          }
          
          # remove the contribution to the compartment counts in the configuration matrix
          .rows_to_update <- 1:epimodel$.ind_final_config
          
          epimodel$config_mat[, epimodel$states][cbind(.rows_to_update, epimodel$config_mat[.rows_to_update, .subj_ID])]<-epimodel$config_mat[,epimodel$states][cbind(.rows_to_update, epimodel$config_mat[.rows_to_update, .subj_ID])]- 1
          
          # compute the likelihood of the subject's trajectory from a
          # time-inhomogeneous CTMC - only computed when current path is saved
          if(save_path) {
                    
                    # check to see if any additional irms are needed.
                    # if so, check_irm will instatiate the required
                    # matrices and their eigen decompositions
                    check_irm(epimodel)
                    
                    epimodel$likelihoods$subj_likelihood_cur <- calc_subj_likelihood(epimodel = epimodel, subject = subject, subj_ID = .subj_ID, log = TRUE)
          }
          
          # if there were transitions, remove the rows, update the indexing variables
          if(length(.subj_inds) != 0) {
                    
                    # remove the rows relating to the trajectory from the configuration matrix
                    epimodel$config_mat[.subj_inds,] <- NA
                    epimodel$config_mat <- epimodel$config_mat[order(epimodel$config_mat[,"time"]),]
                    
                    # set the index of the final configuration
                    epimodel$.ind_final_config <- epimodel$.ind_final_config - length(.subj_inds)
                   
          }
          
          # update the compartment counts in the observation matrix and the
          # indices for observation times in the configuration matrix
          update_obs_mat(epimodel)
          
}