#' Construct an environment of rate matrices and another environment for the
#' corresponding eigen decompositions.
#' 
#' Each rate matrix in the array has dimension \code{num_states x num_states}.
#' 
#' @inheritParams simulate_epimodel
#'   
#' @return environment with rate matrices
#'   
build_irm <- function(epimodel) {
          
          # if there is no .irm environment within .epimodel, instatiate it and
          # the keys
          if(is.null(epimodel$.irm)) {
                    
                    # get compartments names that will be used to generate keys
                    # and index into the rate matrix environment
                    .unique_configs <- unique(epimodel$config_mat[1:epimodel$.ind_final_config, epimodel$index_states, drop = FALSE])
                    epimodel$.irm_keys <- apply(.unique_configs, 1, paste0, collapse = epimodel$index_states)
                    
                    # instatiate the .irm and .eigen environments
                    epimodel$.irm       <- new.env(hash = TRUE, size = length(epimodel$.irm_keys))
                    epimodel$.eigen     <- new.env(hash = TRUE, size = length(epimodel$.irm_keys))
                    
                    # build rate matrices
                    for(r in 1:length(epimodel$.irm_keys)) {
                              
                              # instatiate the rate matrix for each key
                              epimodel$.irm[[epimodel$.irm_keys[[r]]]] <- build_rate_mat(rates = epimodel$rates, state = .unique_configs[r,, drop = FALSE], params = epimodel$params, flow = epimodel$flow, env = NULL, key = NULL)
                              
                              # instatiate the corresponding eigen decomposition
                              epimodel$.eigen[[epimodel$.irm_keys[[r]]]] <- irm_decomp(irm = epimodel$.irm[[epimodel$.irm_keys[[r]]]], env = NULL, key = NULL)
                              
                    }
                    
          } else { # if an .irm environment exists within .epimodel, merely update the rate matrices in place
                    
                    # get the unique configurations from the keys
                    .unique_configs <- matrix(as.numeric(epimodel$.irm_keys), dimnames = list(NULL,epimodel$index_states))
                    #update the rate matrices and eigen decompositions
                    for(r in 1:length(epimodel$.irm_keys)) {
                              
                              # update the rate matrix
                              build_rate_mat(rates = epimodel$rates, state = .unique_configs[r,, drop = FALSE], params = epimodel$params, flow = epimodel$flow, env = epimodel$.irm, key = epimodel$.irm_keys[[r]]) 
                              
                              # update the eigen decomposition
                              irm_decomp(irm = epimodel$.irm[[epimodel$.irm_keys[[r]]]], env = epimodel$.eigen, key = epimodel$.irm_keys[[r]])
                    }
                    
          }
}