#' Construct an index matrix identifying with which states each rate varies.
#' 
#' Called internally before rates are parsed when an epimodel object is
#' initialized.
#' 
#' @inheritParams init_epimodel
#' 
#' @return binary index matrix indicating which changes in each compartment affect c3

init_rate_map <- function(rates, states) {
          
          # initialize matrix
          rate_map <- matrix(0, nrow = length(rates), ncol = length(states))
          colnames(rate_map) <- states
          
          # identify which states show up as arguments in each rate function
          for(t in 1:length(rates)) {
                    rate_map[t, sapply(states, grepl, rates[t])] <- 1
          }
          
          return(rate_map)
}