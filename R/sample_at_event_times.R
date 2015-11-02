#' Sample the infection status at times of state transition of other subjects
#'
#' @param path vector corresponding to a subject path.
#' @inheritParams draw_trajec
#'
#' @return updated path vector of subject level trajectories
#'
#' @export

sample_at_event_times <- function(path, epimodel) {
        
        # If the model is progressive and has an absorbing state, only need to
        # sample the discrete time skeleton where the endpoints of observation
        # intervals indicate a jump
        if(epimodel$absorbing_states & epimodel$progressive) {
                
                for(t in 1:(epimodel$nobs - 1)) {
                        
                        # draw the status at inter-event times conditional on the
                        # status at observation times only if a change is observed
                        if(path[epimodel$obs_time_inds[t]] != path[epimodel$obs_time_inds[t + 1]]) {
                                
                                # if rates change in the inter-observation time interval,
                                # sample the discrete time skeleton at event times
                                if(diff(epimodel$obs_time_inds[t:(t+1)], lag = 1) != 1) {
                                        
                                        path <- sampleEventSubseq(path, epimodel$tpms, epimodel$tpm_products, epimodel$obs_time_inds[t], epimodel$obs_time_inds[t+1])
                                        
                                }
                                
                        }
                        
                }
                
                #if the model is not progressive but has an absorbing state,
                #no need to deal with the path after absorbtion
        } else if(epimodel$absorbing_states & !epimodel$progressive) {
                
                for(t in 1:(epimodel$nobs - 1)) {
                        
                        # draw the status at inter-event times conditional
                        # on the status at observation times only if the
                        # subject is not in an absorbing state
                        if(!path[epimodel$obs_time_inds[t]] %in% epimodel$absorbing_states) {
                                if(diff(epimodel$obs_time_inds[t:(t+1)], lag = 1) != 1) {
                                        
                                        path <-  sampleEventSubseq(path, epimodel$tpms, epimodel$tpm_products, epimodel$obs_time_inds[t], epimodel$obs_time_inds[t+1])
                                        
                                }
                                
                                
                        } else break
                        
                }
                
                # if the model is both not progressive and lacks an
                # absorbing state, need to sample all inter-event intervals
        } else {
                for(t in 1:(epimodel$nobs - 1)) {
                        
                        if(diff(epimodel$obs_time_inds[t:(t+1)], lag = 1) != 1) {
                                
                                 path <- sampleEventSubseq(path, epimodel$tpms, epimodel$tpm_products, epimodel$obs_time_inds[t], epimodel$obs_time_inds[t+1])
                        }
                        
                }
                
        }
          
          return(path)
}
