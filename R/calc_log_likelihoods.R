#' Calculate likelihood for the population level trajectory and the measurement
#' process.
#' 
#' @inheritParams simulate_epimodel
#' @param pop TRUE/FALSE for whether to compute the population-level
#'   log-likelihood
#' @param obs TRUE/FALSE for whether to compute the likelihood of the
#'   measurement process.
#' @param epimodel_envir TRUE/FALSE for whether the function is called
#'   internally within the epimodel environment; if not, the function explicitly
#'   returns a list with the log-likelihoods that were requested.
#'   
#' @return vector with log-likelihood of the population level trajectory and
#'   measurement process log-likelihood.
#' @export
#' 
calc_log_likelihoods <- function(epimodel, pop = TRUE, obs = TRUE, epimodel_envir = FALSE) {
          
          if(pop == TRUE) {
                    
                    if(epimodel_envir == FALSE && is.null(epimodel$.ind_final_config)) {
                              epimodel$.ind_final_config <- which(epimodel$config_mat[,"time"] == max(epimodel$obstimes))
                    }
                    
                    # compute all rates, event specific rates, and hazards for
                    # the observed events
                    
                    # calculate the rate for each type of transition on the subject level
                    .all_rates          <- do.call(cbind, sapply(epimodel$rates, do.call, list(state = epimodel$config_mat[1:(epimodel$.ind_final_config - 1),], params = epimodel$params)))
                    
                    # extract the appropriate subject level rate for each transition event
                    .event_rates        <- .all_rates[cbind(1:(epimodel$.ind_final_config - 1), epimodel$config_mat[2:epimodel$.ind_final_config, "Event"])]
                    
                    # compute hazards by multiplying the sybject level rates by
                    # the number of subjects at risk for each transition
                    .hazards            <- rowSums(.all_rates * epimodel$config_mat[1:(epimodel$.ind_final_config - 1), epimodel$rate_map])
                    
                    # compute the time intervals
                    .time_diffs         <- diff(epimodel$config_mat[1:epimodel$.ind_final_config, "time"], lag = 1)
                    
                    # sum the log-likelihood for the initial configuration,
                    # exponential waiting time densities, and tail probability
                    pop_log_likelihood <- sum(vapply(X = epimodel$config_mat[1,epimodel$.config_inds], FUN = epimodel$d_initdist, params = epimodel$params, log = TRUE, FUN.VALUE = numeric(1))) + sum(log(.event_rates)) - sum(.hazards * .time_diffs)
                    
          }
          
          if(obs == TRUE) {
                    
                    obs_log_likelihood <- sum(epimodel$d_meas_process(state = epimodel$obs_mat, meas_vars = epimodel$meas_vars, params = epimodel$params, log = TRUE))
                    
          }
          
          # Either modify the log-likelihood objects in the epimodel environment or return a list (default)
          if(epimodel_envir == TRUE) {
                    
                    epimodel$pop_log_likelihood <- pop_log_likelihood
                    epimodel$obs_log_likelihood <- obs_log_likelihood
                    
          } else if(epimodel_envir == FALSE) {
                    
                    return(list(pop_log_likelihood = pop_log_likelihood, obs_log_likelihood = obs_log_likelihood))
                    
          }
}