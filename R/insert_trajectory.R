#' Insert a newly sampled subject-level trajectory so its contribution is 
#' reflected in the configuration matrix.
#' 
#' This function also updates several bookkeeping objects - .ind_final_config
#' 
#' @inheritParams sample_path
#' @param reinsertion TRUE/FALSE for whether the trajectory is to be be
#'   reinserted following a MH rejection
#'   
#' @return updated configuration matrix in the epimodel environment
#' @export

insert_trajectory <- function(epimodel, subject, subj_ID, reinsertion) {
          
          # first re-order the configuration matrix if needed, then copy the
          # preceding rows for the new trajectory
          if(epimodel$.ind_final_config != (epimodel$.subj_row_ind - 1)) {
                    
                    # reorder the matrix
                    epimodel$config_mat[1:(epimodel$.subj_row_ind - 1),] <- epimodel$config_mat[order(epimodel$config_mat[1:(epimodel$.subj_row_ind - 1), "time"]), ]
                              
                    # get the subject indices
                    .subj_inds <- which(epimodel$config_mat[,"ID"] == subject)
                    
                    for(k in seq_along(.subj_inds)) {
                              
                              # # copy the states from the previous event times
                              epimodel$config_mat[.subj_inds[k], epimodel$states] <- epimodel$config_mat[.subj_inds[k] - 1, epimodel$states]
                              
                              # copy the configurations from the previous event times
                              epimodel$config_mat[.subj_inds[k], epimodel$.config_inds[-subject]] <- epimodel$config_mat[.subj_inds[k] - 1, epimodel$.config_inds[-subject]] 
                              
                    }
                    
                    # update index for the final configuration
                    epimodel$.ind_final_config <- epimodel$.subj_row_ind - 1
                    
          }
          
          # compute the likelihood of the subject's trajectory from a 
          # time-inhomogeneous CTMC - only when it is not a reinsertion
          # following a MH rejection
          if(!reinsertion) {
                    
                    # check to see if any additional irms are needed.
                    # if so, check_irm will instatiate the required
                    # matrices and their eigen decompositions
                    check_irm(epimodel)
                    
                    epimodel$likelihoods$subj_likelihood_new <- calc_subj_likelihood(epimodel = epimodel, subject = subject, subj_ID = subj_ID, log = TRUE)
                    
          }
           
          
          # add the subject's contribution back into the compartment counts
          .rows_to_update <- 1:epimodel$.ind_final_config
          
          epimodel$config_mat[, epimodel$states][cbind(.rows_to_update, epimodel$config_mat[.rows_to_update, subj_ID])] <- epimodel$config_mat[,epimodel$states][cbind(.rows_to_update, epimodel$config_mat[.rows_to_update, subj_ID])] + 1
          
          # compute the likelihood for the new population-level trajectory -
          # only when not a reinsertion following a MH rejection
          if(!reinsertion) {
                    
                    epimodel$likelihoods$pop_likelihood_new <- calc_pop_likelihood(epimodel, log = TRUE)
                    
          }
}