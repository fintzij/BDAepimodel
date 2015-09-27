#' Updates the matrix of emission probabilities
#'
#' @param epimodel 
#'
#' @return updated matrix of emission probabilities

build_emission_mat <- function(epimodel) {
          
          epimodel$.emission_mat[-epimodel$unmeasured_vars,]
          
}