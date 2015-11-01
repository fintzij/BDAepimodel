#' Insert a newly sampled subject-level trajectory so its contribution is 
#' reflected in the configuration matrix.
#' 
#' This function also updates several bookkeeping objects - .ind_final_config
#' 
#' @inheritParams sample_path
#' @param reinsertion TRUE/FALSE for whether the trajectory is to be be
#'   reinserted following a MH rejection
#'   
#' @return updated epimodel object
#' @export

insert_trajectory <- function(epimodel, subject, reinsertion) {
          
          # first re-order the configuration matrix if needed, then copy the
          # preceding rows for the new trajectory
          if(epimodel$ind_final_config != (epimodel$subj_row_ind - 1)) {
                    
                    # reorder the matrix
                    row_ord <- order(epimodel$pop_mat[, "time"])
                    
                    epimodel$pop_mat <- reorderMat(epimodel$pop_mat, row_ord) # reorder matrix
                    colnames(epimodel$pop_mat) <- c("time", "ID", "Event", epimodel$states) # replace column names
                    
                    epimodel$config_mat <- reorderMat(epimodel$config_mat, row_ord) # reorder config_mat
                              
                    # get the subject indices
                    subj_inds <- which(epimodel$pop_mat[,"ID"] == subject)
                    
                    for(k in seq_along(subj_inds)) {
                              
                              # copy the states from the previous event times
                              epimodel$pop_mat[subj_inds[k], epimodel$states] <- epimodel$pop_mat[subj_inds[k] - 1, epimodel$states]
                              
                              # copy the configurations from the previous event times
                              epimodel$config_mat[subj_inds[k], -subject] <- epimodel$config_mat[subj_inds[k] - 1, -subject] 
                              
                    }
                    
                    # update index for the final configuration and obs time indices
                    epimodel$ind_final_config <- epimodel$subj_row_ind - 1
                    epimodel$obs_time_inds <- getObsTimeInds(epimodel$pop_mat, epimodel$obstimes)
          }
          
          # compute the likelihood of the subject's trajectory from a 
          # time-inhomogeneous CTMC - only when it is not a reinsertion
          # following a MH rejection
          if(!reinsertion) {
                    
                    # get the vector of keys
                    epimodel$keys <- retrieveKeys(1:epimodel$ind_final_config, epimodel$irm_key_lookup, epimodel$pop_mat, epimodel$index_state_num)
                    
                    epimodel$likelihoods$subj_likelihood_new <- subjectLikelihood(subject, pop_mat = epimodel$pop_mat, config_mat = epimodel$config_mat, irm_array = epimodel$irm, initdist = epimodel$initdist, keys = epimodel$keys, inds = c(1, c(1:epimodel$ind_final_config)[-epimodel$obs_time_inds], epimodel$ind_final_config), loglik = TRUE)
                    
          }
          
          # add the subject's contribution back into the compartment counts
          resolveSubjContrib(epimodel$pop_mat, 1:epimodel$ind_final_config, epimodel$config_mat[,subject], insertion =  TRUE)
          
          # compute the likelihood for the new population-level trajectory -
          # only when not a reinsertion following a MH rejection
          if(!reinsertion) {
                    
                    # update the keys to reflect the insertion - no dimension
                    # change from previous retrieval, so only update keys where
                    # subject was in an index state
                    index_contrib <- which(epimodel$config_mat[,subject] %in% epimodel$index_state_num)
                    epimodel$keys[index_contrib]<- retrieveKeys(index_contrib, epimodel$irm_key_lookup, epimodel$pop_mat, epimodel$index_state_num)
                    
                    # deprecated - only need to retrieve keys where subject was in an index state
                    # epimodel$keys<- retrieveKeys(1:epimodel$ind_final_config, epimodel$irm_key_lookup, epimodel$pop_mat, epimodel$index_state_num)
                    
                    # check to see if any additional irms are needed. check_irm will
                    # instatiate the required matrices and their eigen decompositions
                    if(any(epimodel$keys == 0)) {
                              
                              epimodel <- append_missing(epimodel)
                              
                    }
                    
                    epimodel$likelihoods$pop_likelihood_new <- populationLikelihood(pop_mat = epimodel$pop_mat, irm_array = epimodel$irm, initdist = epimodel$initdist, initdist_param_inds = epimodel$initdist_param_inds, flow_inds = epimodel$flow_inds, keys = epimodel$keys, inds = c(1, c(1:epimodel$ind_final_config)[-epimodel$obs_time_inds], epimodel$ind_final_config), loglik = TRUE)
                    
          } 
          
          return(epimodel)
}