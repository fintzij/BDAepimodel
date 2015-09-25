#' Calculate likelihood for the population level trajectory and the measurement
#' process.
#' 
#' @inheritParams simulate_epimodel
#' @param pop TRUE/FALSE for whether to compute the population-level
#'   log-likelihood
#' @param obs TRUE/FALSE for whether to compute the likelihood of the
#'   measurement process.
#' @param log TRUE/FALSE for whether to return the log-likelihood or likelihood. 
#'   
#' @return vector with log-likelihood of the population level trajectory and
#'   measurement process log-likelihood.
#' @export
#' 
calc_likelihoods <- function(epimodel, pop = TRUE, obs = TRUE, epimodel_envir = FALSE, log = FALSE) {
          
          if(pop == TRUE) {
                    
                    if(epimodel_envir == FALSE && is.null(epimodel$.ind_final_config)) {
                              epimodel$.ind_final_config <- which(epimodel$config_mat[,"time"] == max(epimodel$obstimes))
                    }
                    
                    # get the appropriate indices to exclude the intermediate observation time rows
                    .left_endpoints <- c(1, which(epimodel$config_mat[,"ID"] != 0))
                    .right_endpoints <- c(.left_endpoints[-1], epimodel$.ind_final_config)
                    
                    # compute all rates, event specific rates, and hazards for
                    # the observed events
                    
                    # calculate the rate for each type of transition on the subject level
                    .all_rates          <- do.call(cbind, sapply(epimodel$rates, do.call, list(state = epimodel$config_mat[.left_endpoints,], params = epimodel$params)))
                    
                    # extract the appropriate subject level rate for each transition event
                    .event_rates        <- .all_rates[cbind(1:(length(.left_endpoints)-1), epimodel$config_mat[.right_endpoints[-length(.right_endpoints)], "Event"])]
                    
                    # compute hazards by multiplying the sybject level rates by 
                    # the number of subjects at risk for each transition -
                    # indicated by the state_lookup vector
                    .hazards            <- rowSums(.all_rates * epimodel$config_mat[.left_endpoints, epimodel$state_lookup])
                    
                    # compute the time intervals
                    .time_diffs         <- diff(epimodel$config_mat[c(.left_endpoints, epimodel$.ind_final_config), "time"], lag = 1)
                    
                    # sum the log-likelihood for the initial configuration,
                    # exponential waiting time densities, and tail probability
                    pop_likelihood <- sum(vapply(X = epimodel$config_mat[1,epimodel$.config_inds], FUN = epimodel$d_initdist, FUN.VALUE = numeric(1), params = epimodel$params, log = TRUE)) + sum(log(.event_rates)) - sum(.hazards * .time_diffs)
                    
                    if(log == FALSE) {
                              pop_likelihood <- exp(pop_likelihood)
                    }
                    
          }
          
          if(obs == TRUE) {
                    
                    obs_likelihood <- sum(epimodel$d_meas_process(state = epimodel$obs_mat, meas_vars = epimodel$meas_vars, params = epimodel$params, log = TRUE))
                    
                    if(log == FALSE) {
                              obs_likelihood <- exp(obs_likelihood)
                    }
                    
          }
          
          # Either modify the log-likelihood objects in the epimodel environment or return a list (default)
          if(epimodel_envir == TRUE) {
                    
                    if(pop & obs) {
                              epimodel$pop_likelihood <- pop_likelihood
                              epimodel$obs_likelihood <- obs_likelihood
                              
                    } else if(pop & !obs) {
                              epimodel$pop_likelihood <- pop_likelihood
                              
                    } else if(!pop & obs) {
                              epimodel$obs_likelihood <- obs_likelihood
                              
                    }
                    
          } else if(epimodel_envir == FALSE) {
                    
                    if(pop & obs) {
                              return(list(pop_likelihood = pop_likelihood, obs_likelihood = obs_likelihood))
                              
                    } else if(pop & !obs) {
                              return(pop_likelihood)
                              
                    } else if(!pop & obs) {
                              return(obs_likelihood)
                              
                    }
                    
          }
}