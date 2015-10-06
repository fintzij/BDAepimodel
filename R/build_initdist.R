#' Construct the vector with initial state probabilities for a single subject at
#' time t0.
#' 
#' @inheritParams simulate_epimodel
#'   
#' @return vector of initial state probabilities for a single subject
#' @export

build_initdist <- function(epimodel) {
          
          initdist <- epimodel$params[epimodel$initdist_params]
          
          return(initdist)
}