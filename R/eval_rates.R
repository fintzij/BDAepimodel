#' Evaluate the rate functions at the current state and the parameter values,
#' returning a matrix of rates for use in the \code{simulate_epimodel} function.
#' 
#' @param rate_fcns list of rate functions created by extract_rate_fcns
#' @param state current configuration of the population state vector
#' @inheritParams init_epimodel
#'   
#' @return matrix of rates
#'   
get_rate_mat <- function(rate_fcns, state, flow, params) {
          
          rate_mat <- matrix(0, nrow = length(state), ncol = nrow(flow))

          for(k in 1:nrow(flow)){
                    which_state <- which(flow[k,] == -1)
          }
          
}