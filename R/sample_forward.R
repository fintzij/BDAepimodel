#' Forward path sampling.
#' 
#' @param path matrix for storing path, possibly with first transition completed
#' @inheritParams sample_path_in_interval
#' @param init_time
#' @param final_time
#' @param init_state
#' @param final_state
#' @param irm_key
#'   
#' @return updated bookkeeping matrix with a valid path (no update if constant
#'   path)
#'   
#' @export

sample_forward <- function(path, epimodel, subject, init_time, final_time, init_state, final_state, irm_key) {
          
          valid_path          <- FALSE
          ind_init            <- if(is.na(path[1, 1])) {1} else {2}
          
          while(valid_path == FALSE) {
                    
                    keep_going         <- TRUE
                    .t                 <- init_time
                    cur_state          <- init_state
                    ind_cur            <- ind_init
                    
                    while(keep_going == TRUE) {
                              
                              if(cur_state %in% epimodel$absorbing_states) {
                                        
                                        # stop forward sampling
                                        keep_going = FALSE

                                        # determine if the path is valid 
                                        if(cur_state == final_state) {
                                                  
                                                  valid_path <- TRUE
                                        
                                        } else {
                                                  
                                                  valid_path <- FALSE
                                                  path[ind_init : ind_cur, ] <- NA
                                                  
                                                  }
                                                                      
                                        break
                              } 
                              
                              # sample the next transition time
                              .dt <- if(epimodel$irm[cur_state, cur_state, irm_key] == 0) {
                                        Inf
                              } else {
                                        rexp(1, rate = abs(epimodel$irm[cur_state, cur_state, irm_key]))
                              }
                              
                              .t <- .t + .dt
                              
                              # if the next transition time is after the right endpoint, 
                              # stop sampling forward and determine if the path is valid
                              if(.t > final_time) {
                                        
                                        # stop forward sampling
                                        keep_going = FALSE
                                        
                                        # determine if the path is valid. If valid, keep 
                                        # the extra rows and set epimodel$subj_row_ind to 
                                        # ind_cur. otherwise, ind_cur and cur_state will be 
                                        # reset in the next attempt
                                        if(cur_state == final_state) {
                                                  
                                                  valid_path <- TRUE
                                                  
                                        } else {
                                                  valid_path <- FALSE
                                                  path[ind_init : ind_cur, ] <- NA
                                        }
                                        
                              # if .t <= final_time, sample the next state 
                              } else {
                                        
                                        if(ind_cur > nrow(path)) {
                                                  path <- insert_row(path, nrow(path), NA)
                                        }
                                        
                                        # sample the next state
                                        next_state <- .Internal(sample(epimodel$num_states, 1, FALSE, prob = pmax(epimodel$irm[cur_state, , irm_key], 0)))
                                        # next_state <- sample.int(epimodel$num_states, 1, prob = pmax(epimodel$irm[cur_state, , irm_key], 0))
                                        
                                        # identify the event
                                        event <- which((epimodel$flow[, cur_state] == -1) & (epimodel$flow[, next_state] == 1))
                                        
                                        # update the configuration matrix
                                        path[ind_cur, ] <- c(.t, event, next_state)
                                        
                                        # update the sampling object
                                        ind_cur            <- ind_cur + 1
                                        cur_state          <- next_state 
                              }
                    }
          }
          
          if(ind_cur == 1) {
                    return(NULL)
                    
          } else {
                    return(path[1:(ind_cur - 1), , drop = FALSE])
          }
}