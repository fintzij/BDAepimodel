#' Update the list of tpms and tpm products to reflect the removal or insertion 
#' of a subject-level trajectory.
#' 
#' Gets called before \code{remove_trajectory} in the removal direction, and
#' after \code{draw_trajec} in the insertion direction.
#' 
#' @inheritParams draw_trajec
#' @param direction "removal" or "insertion"
#'   
#' @return updated lists of tpms within the .epimodel environment
#' 
#' @export

update_tpms <- function(epimodel, subject, direction) {
          
          # get indices for subject-related events in the configuration matrix
          .subj_inds <- which(epimodel$config_mat[,"ID"] == subject)
          
          if(direction == "removal") {
          
                    # re-order the .tpms and .tpm_products objects
                    .tpm_order <- c(setdiff(1:length(epimodel$.tpms), .subj_inds), .subj_inds)
                    
                    # set the offending matrices to NULL
                    epimodel$.tpms[.subj_inds] <- epimodel$.tpm_products[.subj_inds] <- list(NULL)
                    
                    # re-order the lists to suppress the NULL matrices
                    epimodel$.tpms <- epimodel$.tpms[.tpm_order]
                    epimodel$.tpm_products <- epimodel$.tpm_products[.tpm_order]
                    
                    # Now rebuild the tpm lists in the places that need to be modified
                    # get indices of tpms that need to be rebuilt
                    epimodel$.tpms_to_build <- get_tpms_to_build(epimodel, subject = subject)
                    
                    # update the transition probability matrix sequences
                    build_tpm_seqs(epimodel)
          
                    
          } else if(direction == "insertion") {
                    
                    # re-order the .tpms and .tpm_products objects
                    .tpm_order <- order(c(setdiff(1:length(epimodel$.tpms), .subj_inds), .subj_inds))
                    
                    epimodel$.tpms <- epimodel$.tpms[.tpm_order]
                    epimodel$.tpm_products <- epimodel$.tpm_products[.tpm_order]
                    
                    # Now rebuild the tpm lists in the places that need to be modified
                    # get indices of tpms that need to be rebuilt
                    epimodel$.tpms_to_build <- get_tpms_to_build(epimodel, subject = subject)
                    
                    # update the transition probability matrix sequences
                    build_tpm_seqs(epimodel)
          
          }
}