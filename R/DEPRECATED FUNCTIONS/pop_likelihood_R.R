calc_pop_likelihood <- function(epimodel, params = NULL, log = TRUE) {
          
          if(is.null(epimodel$ind_final_config)) {
                    epimodel$ind_final_config <- which(epimodel$pop_mat[,"time"] == max(epimodel$obstimes))
          }
          
          # get the appropriate indices to exclude the intermediate observation time rows
          left_endpoints     <- c(1, which(epimodel$pop_mat[,"ID"] != 0))
          right_endpoints    <- c(left_endpoints[-1], epimodel$ind_final_config)
          
          # compute all rates, event specific rates, and hazards for
          # the observed events
          
          # calculate the rate for each type of transition on the subject level
          if(is.null(params)) {
                    
                    all_rates <- do.call(cbind, lapply(epimodel$rates, do.call, list(state = epimodel$pop_mat[left_endpoints, epimodel$states], params = epimodel$params)))  
                    
          } else {
                    
                    all_rates<- do.call(cbind, lapply(epimodel$rates, do.call, list(state = epimodel$pop_mat[left_endpoints, epimodel$states], params = params)))
                    
          }
          
          # extract the appropriate subject level rate for each transition event
          event_rates        <- all_rates[cbind(1:(length(left_endpoints)-1), epimodel$pop_mat[right_endpoints[-length(right_endpoints)], "Event"])]
          
          # compute hazards by multiplying the sybject level rates by 
          # the number of subjects at risk for each transition -
          # indicated by the state_lookup vector
          hazards            <- rowSums(all_rates * epimodel$pop_mat[left_endpoints, epimodel$state_lookup])
          
          # compute the time intervals
          time_diffs         <- diff(epimodel$pop_mat[c(left_endpoints, epimodel$ind_final_config), "time"], lag = 1)
          
          # sum the log-likelihood for the initial configuration,
          # exponential waiting time densities, and tail probability
          pop_likelihood <- sum(epimodel$pop_mat[1, epimodel$states[epimodel$initdist_param_inds]] * log(epimodel$initdist[epimodel$initdist_param_inds])) + sum(log(event_rates)) - sum(hazards * time_diffs)
          
          if(log == FALSE) {
                    pop_likelihood <- exp(pop_likelihood)
          }
          
          return(pop_likelihood)
}