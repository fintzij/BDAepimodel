#' Construct a single rate matrix from a vector of rates and a flow matrix.
#' 
#' @param rates vector containing values of evaluated rates
#' @inheritParams init_epimodel
#' @param env environment into which to instatiate the rate matrix
#' @param key the key for the hash index
#' @param update TRUE/FALSE for whether to update an existing rate matrix or to
#'   instatiate a new one
#'   
#' @return rate matrix of dimension \code{num_states x num_states}
#' @export

build_rate_mat <- function(rates, state, params, flow, env = NULL, key = NULL, update = TRUE) {
          
          .rates <- sapply(rates, do.call, list(state = state, params = params))
          
          # if no key is provided, a rate matrix must be instatiated, otherwise
          # modifications occur in-place
          if(is.null(key) || update == FALSE) {
                    
                    .rate_mat <- matrix(0, nrow = ncol(flow), ncol = ncol(flow))
                    
                    for(k in 1:nrow(flow)) {
                              .rate_mat[flow[k,] == -1, flow[k,] == 1] <- .rate_mat[flow[k,] == -1, flow[k,] == 1] + .rates[k]
                              .rate_mat[flow[k,] == -1, flow[k,] == -1] <- .rate_mat[flow[k,] == -1, flow[k,] == -1] - .rates[k]
                    }
                    
                    return(.rate_mat)
                    
          } else {
                    env[[key]][1:ncol(flow), 1:ncol(flow)] <- 0
                    
                    for(k in 1:nrow(flow)) {
                              env[[key]][flow[k,] == -1, flow[k,] == 1] <- env[[key]][flow[k,] == -1, flow[k,] == 1] + .rates[k]
                              env[[key]][flow[k,] == -1, flow[k,] == -1] <- env[[key]][flow[k,] == -1, flow[k,] == -1] - .rates[k]
                    }
          }
          
}