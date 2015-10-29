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
#' @export
#'

expand_config_mat <- function(epimodel, buffer_size = NULL) {

          # insert additional rows for configurations at observation times
          if(!all(epimodel$obstimes %in% epimodel$pop_mat[,"time"])) {

                    # get obstimes to insert
                    inds_to_insert <- which(!epimodel$obstimes %in% epimodel$pop_mat[,"time"])
                    obstimes_to_insert <- epimodel$obstimes[inds_to_insert]

                    for(k in 1:length(obstimes_to_insert)) {

                              insert_ind <- sum(epimodel$pop_mat[,"time"] <= obstimes_to_insert[k]) + 1

                              epimodel$pop_mat <- insert_row(epimodel$pop_mat, insert_ind, epimodel$pop_mat[insert_ind - 1, ])
                              epimodel$pop_mat[insert_ind, c("time", "ID", "Event")] <- c(obstimes_to_insert[k], 0, 0)

                              epimodel$config_mat <- insert_row(epimodel$config_mat, insert_ind, epimodel$config_mat[insert_ind - 1, ])

                    }
          }

          # get the index for the final configuration
          epimodel$ind_final_config <- which(epimodel$pop_mat[,"time"] == max(epimodel$obstimes))

          # set the vector of observation time indices
          epimodel$obs_time_inds <- which(epimodel$pop_mat[,"ID"] == 0)

          # if the column of irm_keys has not been generated, generate it

          # add buffer of NAs
          if(is.null(buffer_size)) {
                    buffer_size <- floor(0.3 * epimodel$ind_final_config)
          }

          # expand the configuration matrix
          epimodel$pop_mat <- rbind(epimodel$pop_mat, matrix(NA, nrow = buffer_size, ncol = ncol(epimodel$pop_mat)))
          epimodel$config_mat <- rbind(epimodel$config_mat, matrix(NA, nrow = buffer_size, ncol = ncol(epimodel$config_mat)))
          
          if(!is.null(epimodel$tpms)) {
                    
                    # create empty tpm matrices
                    new_tpm_array  <- array(0, dim = c(length(epimodel$states), length(epimodel$states), buffer_size))
                    new_prod_array <- array(0, dim = c(length(epimodel$states), length(epimodel$states), buffer_size))
                    
                    # append to existing tpm array
                    epimodel$tpms <- joinCubes(epimodel$tpms, new_tpm_array)
                    epimodel$tpm_products <- joinCubes(epimodel$tpm_products, new_prod_array)
          }

          return(epimodel)
}