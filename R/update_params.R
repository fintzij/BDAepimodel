#' Replace the vector of parameters in an epimodel object with new values and 
#' optionally update the population level likelihood or the measurement process 
#' likelihood.
#' 
#' @param epimodel
#' @param params
#' @param pop_likelihood
#' @param obs_likelihood
#'   
#' @return updated parameter vector and likelihood objects in the epimodel
#'   environment
#'   
#' @export

update_params <- function(epimodel, params, pop_likelihood = NULL, obs_likelihood = NULL) {
          
          epimodel$params <- params
          
          if(!is.null(pop_likelihood)) {
                    epimodel$likelihoods$pop_likelihood_cur <- pop_likelihood
                    
          }
          
          if(!is.null(obs_likelihood)) {
                    epimodel$likelihoods$obs_likelihood <- obs_likelihood
          }
}