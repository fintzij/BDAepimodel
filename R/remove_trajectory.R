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

          # get indices related to the subject
          subj_inds          <- which(epimodel$pop_mat[,"ID"] == subject)
          
          # save the current path if not removing it b/c of a M-H rejection
          if(save_path) {
                    epimodel$path_cur  <- cbind(epimodel$pop_mat[c(1, subj_inds, epimodel$ind_final_config), c("time", "Event")], epimodel$config_mat[c(1, subj_inds, epimodel$ind_final_config), subject])
          }
          
          # remove the contribution to the compartment counts in the configuration matrix
          resolveSubjContrib(epimodel$pop_mat, 1:epimodel$ind_final_config, epimodel$config_mat[,subject], insertion =  FALSE)
          
          # compute the likelihood of the subject's trajectory from a
          # time-inhomogeneous CTMC - only computed when current path is saved
          if(save_path) {
                    
                    # update keys
                    epimodel$keys <- retrieveKeys(1:epimodel$ind_final_config, epimodel$irm_key_lookup, epimodel$pop_mat, epimodel$index_state_num)
          
                    # check to see if any additional irms are needed. check_irm will
                    # instatiate the required matrices and their eigen decompositions
                    if(any(epimodel$keys == 0)) {
                    
                        epimodel <- append_missing(epimodel)
                    
                    }
          
                    epimodel$likelihoods$subj_likelihood_cur <- subjectLikelihood(subject, pop_mat = epimodel$pop_mat, config_mat = epimodel$config_mat, irm_array = epimodel$irm, initdist = epimodel$initdist, keys = epimodel$keys, inds = c(1, c(1:epimodel$ind_final_config)[-epimodel$obs_time_inds], epimodel$ind_final_config), loglik = TRUE)
          }

          # if there were transitions, remove the rows, update the indexing variables
          if(length(subj_inds)) {
          
                  # remove the rows relating to the trajectory from the configuration matrix
                  epimodel$pop_mat[subj_inds, ] <- NA # set offending rows to NA
                  epimodel$keys <- epimodel$keys[-subj_inds, , drop = FALSE]
                  
                  row_ord <- order(epimodel$pop_mat[,"time"]) # get new order
                  
                  epimodel$pop_mat <- reorderMat(epimodel$pop_mat, row_ord) # reorder matrix
                  colnames(epimodel$pop_mat) <- c("time", "ID", "Event", epimodel$states) # replace column names
                  
                  epimodel$config_mat[subj_inds,] <- NA # do the same with the configuration matrix
                  epimodel$config_mat <- reorderMat(epimodel$config_mat, row_ord)
                  
                  # set the index of the final configuration
                  epimodel$ind_final_config <- epimodel$ind_final_config - length(subj_inds)
          
          }

          # update the compartment counts in the observation matrix and the
          # indices for observation times in the configuration matrix
          epimodel$obs_time_inds <- getObsTimeInds(epimodel$pop_mat, epimodel$obstimes)
          epimodel$obs_mat[, epimodel$meas_vars_aug] <- epimodel$pop_mat[epimodel$obs_time_inds, epimodel$meas_vars]
          
          return(epimodel)
}