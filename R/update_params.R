#' Replace the vector of parameters in an epimodel object with new values and 
#' update the population level likelihood, the measurement process 
#' likelihood, and the eigen decomposition objects for irms. 
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
          
          return(epimodel)
}