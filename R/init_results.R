#' Initialize list for storing results
#'
#' @inheritParams simulate_epimodel
#' @param sim_settings bookkeeping list of simulation settings
#'
#' @return list containing the matrix of posterior parameter draws, configurations, log-likelihoods, total simulation time, and seed. 
#'
init_results <- function(epimodel) {
          
          sim_settings <- epimodel$sim_settings
          
          # vector for storing log-likelihood
          log_likelihood <- rep(0, floor(sim_settings$niter / sim_settings$params_every))
          
          # matrix for storing parameter updates
          params <- matrix(0, nrow = floor(sim_settings$niter / sim_settings$params_every), ncol = length(epimodel$params))
          colnames(params) <- names(epimodel$params)
          
          # list for storing configurations
          configs <- vector(mode = "list", length = floor(sim_settings$niter / sim_settings$configs_every))
          
          results <- list(time = NULL, seed = NULL, log_likelihood = log_likelihood, params = params, configs = configs)
}