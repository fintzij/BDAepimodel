#' Insert a lookup table to connect state indicators to state labels. 
#'
#' @inheritParams simulate_epimodel 
#'
#' @return data.frame for lookup table for state names and state codes.

init_state_lookup <- function(states) {
          
          state_lookup <- data.frame(state_name = states,
                                              state_code = 1:length(states))
          
          return(state_lookup)
}