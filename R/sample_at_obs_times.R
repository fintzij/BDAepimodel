#' Sample the infection status for the specified subject at observation times.
#'
#' @inheritParams draw_trajec
#'
#' @return updated configuration matrix
#'
#' @export

sample_at_obs_times <- function(config_mat, epimodel, subject) {
        
        # initialize the vector
        state_vec <- epimodel$obstimes
        
        # draw the last state
        state_vec[epimodel$nobs] <- prev_state <- sample.int(n = epimodel$num_states, size = 1, prob = colSums(epimodel$fb_mats[,,epimodel$nobs - 1]))
        
        # backward pass to draw each previous state in succession
        for(s in (epimodel$nobs - 1) : 1) {
                state_vec[s] <- prev_state <- sample.int(n = epimodel$num_states, size = 1, prob = epimodel$fb_mats[, prev_state, s])
                
        }
        
        # update the configuration matrix
        config_mat[epimodel$obs_time_inds, subj_ID] <- state_vec
        
        # if the model is absorbing and progressive, fix the state in constant
        # inter-observation intervals and set the state after absorbtion
        if(epimodel$progressive & epimodel$absorbing_states) {
                
                # identify the first observation time when the subject is in
                # the absorbing state, possibly never
                if(any(state_vec %in% epimodel$absorbing_states)) {
                        absorb_ind <- which(state_vec %in% epimodel$absorbing_states)[1]
                        
                        config_mat[epimodel$obs_time_inds[absorb_ind] : epimodel$ind_final_config, subj_ID] <- state_vec[absorb_ind]
                }
                
                # set the path to be constant in inter-observation intervals
                # with identical endpoints
                for(t in 1:(epimodel$nobs - 1)) {
                        
                        # check if the endpoints of the interval are equal
                        if(state_vec[t] == state_vec[t+1]) {
                                
                                # if so, set the intermediate states to be equal
                                config_mat[epimodel$obs_time_inds[t] : epimodel$obs_time_inds[t+1], subj_ID] <- state_vec[t]
                                
                        }
                        
                }
                
                # if the model merely has an absorbing state but is not progressive,
                # fix the state after absorbtion
        } else if(epimodel$absorbing_states & !epimodel$progressive) {
                
                # identify the first observation time when the subject is in
                # the absorbing state, possibly never
                if(any(state_vec %in% epimodel$absorbing_states)) {
                        absorb_ind <- which(state_vec %in% epimodel$absorbing_states)[1]
                        
                        config_mat[epimodel$obs_time_inds[absorb_ind] : epimodel$ind_final_config, subj_ID] <- state_vec[absorb_ind]
                }
                
                
        }
        
        return(config_mat)
}
