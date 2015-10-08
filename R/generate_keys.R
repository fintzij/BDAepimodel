#' Generate the keys for indexing into the .irm and .eigen environments.
#' 
#' @inheritParams simulate_epimodel
#' @param inds optional vector with row numbers for which keys are desired.
#' @param lookup logical indicating whether to look up the keys in the
#'   .irm_key_lookup object
#'   
#' @return character vector of keys
#' @export

generate_keys <- function(epimodel, inds = NULL, lookup = TRUE) {
          
          # if the indices were not specified, default to
          # 1:(epimodel$.ind_final_config - 1)
          if(is.null(inds)) {
                    inds <- 1:(epimodel$.ind_final_config - 1)
          }
          
          if(lookup) {
                    
                    .key_inds <- apply(epimodel$config_mat[inds, epimodel$index_states, drop = FALSE], 1, match_row, table = epimodel$.irm_key_lookup, return_value = "index")
                    
                    .keys <- rownames(epimodel$.irm_key_lookup[.key_inds, , drop = FALSE])
                    
                    
          } else {
                    .keys <- apply(epimodel$config_mat[inds, epimodel$index_states, drop = FALSE], 1, paste0, collapse = epimodel$index_states)
                    
          }
          
          return(.keys)
}