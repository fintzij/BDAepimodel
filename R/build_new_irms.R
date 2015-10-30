#' Build a new array of rate matrices
#'
#' @param epimodel 
#' @param params names vector of parameters
#'
#' @return epimodel object with updated rate array
#' @export
#'

build_new_irms <- function(epimodel, params) {
          
          # compute new rates
          rates <- do.call(cbind, lapply(epimodel$rates, do.call, list(state = epimodel$irm_key_lookup, params = params)))
          
          # update the rate matrices
          buildRateArray(irm_array = epimodel$irm, rates = rates, flow_inds = epimodel$flow_inds)
          
          return(epimodel)
}