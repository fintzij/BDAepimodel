#' Insert a lookup table to map subjects at risk for each transition to the appropriate rate code label 
#'
#' @inheritParams simulate_epimodel 
#'
#' @return character vector for mapping state codes to compartment names.

init_rate_map <- function(epimodel) {
          
          # matrix to map the state code to the compartment count
          rate_map <- rep(0, nrow(epimodel$flow))
          for(k in 1:length(rate_map)) {
                    rate_map[k] <- colnames(epimodel$flow)[which(epimodel$flow[k,]==-1)]
          }
          
          return(rate_map)
}