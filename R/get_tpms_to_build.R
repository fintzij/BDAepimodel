#' Get the indices for the transition probability matrices to be rebuilt
#'
#' @inheritParams simulate_epimodel
#' @param subject_ID string ID for subject to be removed
#'
#' @return vector of indices
#' @export

get_tpms_to_build <- function(epimodel, subject = NULL) {
          
          # if no subject is specified, all of the tpms must be updated
          if(is.null(subject)) {
                    
                    .tpms_to_build <- 1:(epimodel$.ind_final_config - 1)
                    
          } else { # otherwise, we need to compute which intervals are affected by a subject
                    
                    .subj_ID <- paste0(".X", subject)
                    
                    # first the indices for which the subject was in an index
                    # state that contributes to the rates of other subjects
                    # (e.g. infected in SIR)
                    .inds_index_states <- which(epimodel$config_mat[, .subj_ID][1:(epimodel$.ind_final_config - 1)] %in% epimodel$index_state_num)
                    
                    # then the indices for the intervals just prior to a 
                    # transition for the subject that are not in the index state
                    # indices for that subject. Also get the intervals for the
                    # transitions. 
                    .transition_inds <- c(which(epimodel$config_mat[,"ID"] == subject),
                                          which(diff(epimodel$config_mat[,.subj_ID], lag = 1) != 0))
                    .transition_inds <- .transition_inds[!.transition_inds %in% .inds_index_states]

                    # combine the two vectors
                    .tpms_to_build <- c(.inds_index_states, .transition_inds)
                    
          }
          
          return(.tpms_to_build)
}