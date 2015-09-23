#' Insert a lookup table to map subjects at risk for each transition to the appropriate rate code label 
#'
#' @inheritParams simulate_epimodel 
#'
#' @return character vector for mapping state codes to compartment names.

init_state_lookup <- function(epimodel) {
          
          # matrix to map the state code to the compartment count
          state_lookup <- rep(0, nrow(epimodel$flow))
          for(k in 1:length(state_lookup)) {
                    state_lookup[k] <- colnames(epimodel$flow)[which(epimodel$flow[k,]==-1)]
          }
          
          return(state_lookup)
}