#' Sample the infection status for the specified subject at observation times.
#'
#' @param path vector giving a subject path.
#' @inheritParams draw_trajec
#'
#' @return updated path vector
#'
#' @export

sample_at_obs_times <- function(path, epimodel) {
        
        # initialize the vector
        state_vec <- epimodel$obs_time_inds
        
        # draw the last state
        state_vec[epimodel$nobs] <- prev_state <- .Internal(sample(epimodel$num_states, size = 1, FALSE, colSums(epimodel$fb_mats[,,epimodel$nobs - 1])))
        
#         state_vec[epimodel$nobs] <- prev_state <- sample.int(n = epimodel$num_states, size = 1, prob = colSums(epimodel$fb_mats[,,epimodel$nobs - 1]))
        
        # backward pass to draw each previous state in succession
        for(s in (epimodel$nobs - 1) : 1) {
                  state_vec[s] <- prev_state <- .Internal(sample(epimodel$num_states, 1, FALSE, epimodel$fb_mats[, prev_state, s]))
#                 state_vec[s] <- prev_state <- sample.int(n = epimodel$num_states, size = 1, prob = epimodel$fb_mats[, prev_state, s])
                
        }
        
        # update the path
        path[1:epimodel$ind_final_config] <- state_vec[findInterval(epimodel$pop_mat[1:epimodel$ind_final_config,"time"], epimodel$obstimes)]
        
        return(path)
}
