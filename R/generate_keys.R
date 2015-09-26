#' Generate the keys for indexing into the .irm and .eigen environments. 
#'
#' @inheritParams simulate_epimodel 
#' @param inds optional vector with row numbers for which keys are desired.
#'
#' @return character vector of keys

generate_keys <- function(epimodel, inds = NULL) {
          
          # if the indices were not specified, default to
          # 1:(epimodel$.ind_final_config - 1)
          if(is.null(inds)) {
                    inds <- 1:(epimodel$.ind_final_config - 1)
          }
          
          .keys <- apply(epimodel$config_mat[inds, epimodel$index_states, drop = FALSE], 1, paste0, collapse = epimodel$index_states)
          
          return(.keys)
}