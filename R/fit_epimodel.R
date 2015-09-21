#' Main wrapper to fit a stochastic epimodel using the Bayesian data
#' augmentation algorithm.
#' 
#' @inheritParams simulate_epimodel
#'   
#' @return list containing the epimodel object and a results list. 
#' @export
#' 
fit_epimodel <- function(epimodel) {
          
          # check that the simulation settings have been set
          if(is.null(epimodel$sim_settings)) {
                    stop(sQuote("sim_settings"), "must be specified in the epimodel object.")
          }
          
          if(is.null(epimodel$config_mat)) {
                    warning("An initial configuration was not provided so one has been generated.")
                    epimodel$config_mat
          }
          
          # convert epimodel list to environment to enable multiple assignment
          .epimodel <- list2env(epimodel, parent = emptyenv(), hash = TRUE)
          
          # move simulation settings to internal objects
          .niter              <- .epimodel$sim_settings$niter
          .burnin             <- .epimodel$sim_settings$burnin
          .params_every       <- .epimodel$sim_settings$params_every
          .configs_every      <- .epimodel$sim_settings$configs_every
          .kernel             <- .epimodel$sim_settings$kernel
          .cov_mtx            <- .epimodel$sim_settings$cov_mtx
          .configs_to_redraw  <- .epimodel$sim_settings$configs_to_redraw
          .config_replacement <- ifelse(.configs_to_redraw > .epimodel$popsize, TRUE, FALSE)
          .parallelize        <- .epimodel$sim_settings$parallelize
          .to_estimation_scale <- .epimodel$sim_settings$to_estimation_scale
          .from_estimation_scale <- .epimodel$sim_settings$from_estimation_scale
          
          # initialize list for storing results
          .results <- init_results(.epimodel)
          
          # add buffer to the configuration matrix, also instatiates .epimodel$.ind_final_config     
          .epimodel$config_mat <- expand_config_mat(.epimodel)
          
          # generate .niter parameter samples
          for(k in 1:.niter) {
                    
                    # choose which subjects should be redrawn
                    .subjects <- sample.int(n = .epimodel$popsize, size = .configs_to_redraw, replace = .config_replacement)
                    
                    # re-compute the array of rate matrices
                    .irm <- build_irm(.epimodel)
                    
                    # re-compute the population level likelihood, measurement
                    # process likelihood, and complete data log-likelihood
                    .epimodel$log_likelihoods <- calc_log_likelihoods(epimodel = .epimodel)
                    
                    for(j in 1:.configs_to_redraw) {
                              
                              
                              
                    }
                    
          }
          
          results <- init_results(epimodel)
}