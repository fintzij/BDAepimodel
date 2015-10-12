#' Initialize list for storing results
#' 
#' @inheritParams simulate_epimodel
#' @param sim_settings bookkeeping list of simulation settings
#'   
#' @return list containing the matrix of posterior parameter draws,
#'   configurations, log-likelihoods, total simulation time, and seed.
#' @export
#' 
init_results <- function(epimodel) {
          
          sim_settings <- epimodel$sim_settings
          
          # vector for storing log-likelihood
          log_likelihood <- rep(0, epimodel$sim_settings$niter %/% epimodel$sim_settings$save_params_every)
          
          # matrix for storing parameter updates
          params <- matrix(0, nrow = epimodel$sim_settings$niter %/% epimodel$sim_settings$save_params_every, ncol = length(epimodel$params))
          colnames(params) <- names(epimodel$params)
          
          # matrix for storing acceptance rates
          accepts <- matrix(0, nrow = epimodel$sim_settings$niter - 1, ncol = length(epimodel$params) + 1)
          colnames(accepts) <- c("paths", names(epimodel$params))
          
          # list for storing configurations
          configs <- vector(mode = "list", length = epimodel$sim_settings$niter %/% epimodel$sim_settings$save_configs_every)
          
          results <- list(time = NULL, seed = NULL, log_likelihood = log_likelihood, params = params, accepts = accepts, configs = configs)
}