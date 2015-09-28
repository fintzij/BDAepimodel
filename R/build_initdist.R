#' Construct the vector with initial state probabilities for a single subject at
#' time t0.
#' 
#' @inheritParams simulate_epimodel
#'   
#' @return vector of initial state probabilities for a single subject

build_initdist <- function(epimodel) {
          
          initdist <- vapply(X = c(1,2,3), FUN = epimodel$d_initdist, FUN.VALUE = numeric(1), params = epimodel$params, log = FALSE)
          
          return(initdist)
}