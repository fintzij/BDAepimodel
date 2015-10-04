#' Checks that all of the IRMs that are required exist in the .irm environment 
#' by checking keys against the .irm object list.
#' 
#' @inheritParams simulate_epimodel
#'   
#' @return Instatiated rate matrix in the .irm environment and its eigen
#'   decomposition in the .eigen environment
#' @export

check_irm <- function(epimodel) {
          
          # get the required keys from the configuration matrix
          .unique_configs <- unique(epimodel$config_mat[1:epimodel$.ind_final_config, epimodel$index_states, drop = FALSE])
          .keys <- apply(.unique_configs, 1, paste0, collapse = epimodel$index_states)
          
          # indices of .keys vector for which keys are missing
          .missing_inds <- which(is.na(charmatch(.keys, epimodel$.irm_keys)))
          
          if(length(.missing_inds) != 0) {
                    
                    # missing keys
                    .missing_keys <- .keys[.missing_inds]
                    
                    for(r in 1:length(.missing_keys)) {
                              
                              # add the key to the vector of keys
                              epimodel$.irm_keys <- c(epimodel$.irm_keys, .missing_keys[r])
                              
                              # build the rate matrix
                              epimodel$.irm[[.missing_keys[r]]] <- build_rate_mat(rates = epimodel$rates, state = .unique_configs[.missing_inds[r],, drop = FALSE], params = epimodel$params, flow = epimodel$flow, key = .missing_keys[r], env = epimodel$.irm, update = FALSE)
                              
                              # get the eigen decomposition
                              epimodel$.eigen[[.missing_keys[r]]] <- irm_decomp(irm = epimodel$.irm[[.missing_keys[r]]], key = .missing_keys[r], env = epimodel$.eigen, update = FALSE)
                              
                    }
                    
          }
} 