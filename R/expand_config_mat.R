#' Add buffer rows to a config_mat matrix to enable in-place substitution of
#' configurations.
#' 
#' Each buffer row is made up of NA entries that are inserted at the end of the configuration matrix. This function is designed to operate within an epimodel environment and therefore makes multiple replacements within that environment. 
#' 
#' @inheritParams simulate_epimodel
#' @param buffer_size number of buffer rows to insert, defaulting to 30% of the
#'   number of rows in the current configuration matrix unless otherwise
#'   specified.
#'   
#' @return new configuration matrix and index for the last row in the expanded configuration matrix (corresponding to the configuration at tmax)

expand_config_mat <- function(epimodel, buffer_size = NULL) {
          
          epimodel$.ind_final_config <- which(epimodel$config_mat[,"time"] == max(epimodel$obstimes))
          
          if(is.null(buffer_size)) {
                    buffer_size <- floor(0.3 * epimodel$.ind_final_config)
          }
          
          epimodel$config_mat <- rbind(epimodel$config_mat, matrix(NA, nrow = buffer_size, ncol = ncol(epimodel$config_mat)))
          
}