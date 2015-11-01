#' Draw a path within an interval from an endpoint conditioned homogeneous CTMC.
#' 
#' @inheritParams draw_trajec
#' @param interval index for the interval within which path is to be drawn
#'   
#' @return updated configuration matrix with subject path (excluding at time 0
#'   and tmax) at the end
#'   
#' @export
#' 

sample_path_in_interval <- function(epimodel, subject, interval) {
          
          # set the initial and final states, and times of interval endpoints
          init_state          <- epimodel$config_mat[interval, subject]
          final_state         <- epimodel$config_mat[interval + 1, subject]
          
          init_time           <- epimodel$pop_mat[interval, "time"]
          final_time          <- epimodel$pop_mat[interval + 1, "time"]
          
          irm_key             <- epimodel$keys[interval]
          
          path                <- matrix(nrow = 5, ncol = 3)
          
          # if the endpoints are the same, use the forward sampling
          # algorithm until an appropriate path is generated
          if(init_state == final_state) {
                    
                    path <- sample_forward(path = path,
                                           epimodel = epimodel, 
                                           subject = subject, 
                                           init_time = init_time, 
                                           final_time = final_time, 
                                           init_state = init_state,
                                           final_state = final_state,
                                           irm_key = irm_key)
                    
          } else {
                    # at least one change occurs, so sample it conditionally
                    init_time <- init_time + log(1 - runif(1) * (1 - exp((final_time - init_time) * epimodel$irm[init_state, init_state, irm_key]))) / epimodel$irm[init_state, init_state, irm_key]
                    
                    # choose the next state
                    next_state <- .Internal(sample(epimodel$num_states, 1, FALSE, pmax(epimodel$irm[init_state, , irm_key], 0)))
                    # next_state <- sample.int(epimodel$num_states, 1, prob = pmax(epimodel$irm[init_state, , irm_key], 0))
                    
                    # identify the event and reset the initial state
                    event <- which((epimodel$flow[, init_state] == -1) & (epimodel$flow[, next_state] == 1))
                    init_state <- next_state
                    
                    # insert the first change into the path
                    path[1, ] <- c(init_time, event, init_state)
                    
                    if(!init_state %in% epimodel$absorbing_states) {
                              # sample forward to complete the path
                              path <- sample_forward(path = path,
                                                     epimodel = epimodel, 
                                                     subject = subject,
                                                     init_time = init_time, 
                                                     final_time = final_time,
                                                     init_state = init_state,
                                                     final_state = final_state,
                                                     irm_key = irm_key)         
                    } else {
                              path <- path[1, , drop = FALSE]
                    }
          }
          
          return(path)
}