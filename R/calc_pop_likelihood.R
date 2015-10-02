#' Calculate likelihood for the population level trajectory
#' 
#' @inheritParams simulate_epimodel
#' @param log TRUE/FALSE for whether to return the log-likelihood or likelihood.
#'   
#' @return vector with log-likelihood of the population level trajectory and 
#'   measurement process log-likelihood.
#'   
#' @export
#' 

calc_pop_likelihood <- function(epimodel, log = FALSE) {

          if(is.null(epimodel$.ind_final_config)) {
                    epimodel$.ind_final_config <- which(epimodel$config_mat[,"time"] == max(epimodel$obstimes))
          }
          
          # get the appropriate indices to exclude the intermediate observation time rows
          .left_endpoints     <- c(1, which(epimodel$config_mat[,"ID"] != 0))
          .right_endpoints    <- c(.left_endpoints[-1], epimodel$.ind_final_config)
          
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
          
          return(pop_likelihood)
}