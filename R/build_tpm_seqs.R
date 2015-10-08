#' Construct the sequences of transition probability matrix products.
#' 
#' Updates the .tpms and .tpm_products lists. Resets the .tpms_to_compute vector
#' to be all FALSE. The elements of this vector needing to be updated are then 
#' set to true as each subject level trajectory is re-drawn, and when parameters
#' are updated.
#' 
#' @inheritParams simulate_epimodel
#'   
#' @return modify the .tpms and .tpm_products lists in the .epimodel environment
#'   and reset the .tpms_to_build vector
#' @export  
#' 

build_tpm_seqs <- function(epimodel) {
          
          # get the keys for the tpms that need rebuilding
          .keys <- generate_keys(epimodel, inds = 1:(epimodel$.ind_final_config - 1), lookup = TRUE)
          
          # construct the tpms for all indices in the .tpms_to_build vector
          for(k in seq_along(.keys)) {
                    epimodel$.tpms[[k]] <- build_tpm(values = epimodel$.eigen[[.keys[k]]]$values,
                                                vectors = epimodel$.eigen[[.keys[k]]]$vectors,
                                                inv_vectors = epimodel$.eigen[[.keys[k]]]$inv_vectors,
                                                t0 = epimodel$config_mat[k, "time"],
                                                t1 = epimodel$config_mat[k + 1, "time"])
          }
          
          # get the observation intervals to which the tpms belong
          .intervals <- findInterval(epimodel$config_mat[,"time"], epimodel$obstimes, all.inside = TRUE)
          
          # construct the tpm products in reverse order within each observation interval
          for(k in unique(.intervals)) {
                    
                    # get the indices for all tpms associated with an interval
                    .left_endpoint      <- epimodel$obstimes[k]
                    .right_endpoint     <- epimodel$obstimes[k+1]
                    
                    .subseq_inds <- .Internal(which((epimodel$config_mat[1:(epimodel$.ind_final_config-1),"time"] < .right_endpoint) & (epimodel$config_mat[1:(epimodel$.ind_final_config-1),"time"] >= .left_endpoint)))
                    
                    # Successive multiplication beginning with the rightmost matrix
                    epimodel$.tpm_products[.subseq_inds] <- Reduce("%*%", epimodel$.tpms[.subseq_inds], accumulate = TRUE, right = TRUE)
                    
          }

}