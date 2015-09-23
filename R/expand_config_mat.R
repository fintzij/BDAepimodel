#' Add buffer rows to a config_mat matrix to enable in-place substitution of 
#' configurations.
#' 
#' Each buffer row is made up of NA entries that are inserted at the end of the
#' configuration matrix. This function is designed to operate within an epimodel
#' environment and therefore makes multiple replacements within that
#' environment.
#' 
#' @inheritParams simulate_epimodel
#' @param buffer_size number of buffer rows to insert, defaulting to 30% of the 
#'   number of rows in the current configuration matrix unless otherwise 
#'   specified.
#'   
#' @return new configuration matrix and index for the last row in the expanded
#'   configuration matrix (corresponding to the configuration at tmax)

expand_config_mat <- function(epimodel, buffer_size = NULL) {
          
          # insert additional rows for configurations at observation times
          if(!all(epimodel$obstimes %in% epimodel$config_mat[,"time"])) {
                    
                    # get obstimes to insert
                    .inds_to_insert <- which(!epimodel$obstimes %in% epimodel$config_mat[,"time"])
                    .obstimes_to_insert <- epimodel$obstimes[.inds_to_insert]
                    
                    for(k in 1:length(.obstimes_to_insert)) {
                              
                              .insert_ind <- sum(epimodel$config_mat[,"time"] <= .obstimes_to_insert[k]) + 1
                              
                              epimodel$config_mat <- insert_row(epimodel$config_mat, .insert_ind, epimodel$config_mat[.insert_ind - 1,])
                              epimodel$config_mat[.insert_ind, c("time", "ID", "Event")] <- c(.obstimes_to_insert[k], 0, 0)
                              
                    }
          }
          
          # get the index for the final configuration
          epimodel$.ind_final_config <- which(epimodel$config_mat[,"time"] == max(epimodel$obstimes))
          
          # add buffer of NAs 
          if(is.null(buffer_size)) {
                    buffer_size <- floor(0.3 * epimodel$.ind_final_config)
          }
          
          epimodel$config_mat <- rbind(epimodel$config_mat, matrix(NA, nrow = buffer_size, ncol = ncol(epimodel$config_mat)))
          
}