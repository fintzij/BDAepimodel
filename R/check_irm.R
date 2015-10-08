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
          
          # indices of the configurations which keys are missing
          .missing_inds <- !apply(.unique_configs, 1, match_row, table = epimodel$.irm_key_lookup, return_value = "logical")
          
          if(any(.missing_inds)) {
                    
                    # get configurations for which to generate keys
                    .missing_configs <- .unique_configs[.missing_inds, , drop = FALSE]
                    
                    # generate keys
                    rownames(.missing_configs) <- .missing_keys <- apply(.missing_configs, 1, paste0, collapse = epimodel$index_states)
                    
                    for(r in 1:nrow(.missing_configs)) {
                              
                              # add the key to the vector of keys
                              epimodel$.irm_keys <- c(epimodel$.irm_keys, .missing_keys[r])
                              
                              # build the rate matrix
                              epimodel$.irm[[.missing_keys[r]]] <- build_rate_mat(rates = epimodel$rates, state = .missing_configs[r,, drop = FALSE], params = epimodel$params, flow = epimodel$flow, key = .missing_keys[r], env = epimodel$.irm, update = FALSE)
                              
                              # get the eigen decomposition
                              epimodel$.eigen[[.missing_keys[r]]] <- irm_decomp(irm = epimodel$.irm[[.missing_keys[r]]], key = .missing_keys[r], env = epimodel$.eigen, update = FALSE)
                              
                    }
                    
                    epimodel$.irm_key_lookup <- rbind(epimodel$.irm_key_lookup, .missing_configs)
                    
          }
} 