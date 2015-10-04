#' Forward path sampling.
#' 
#' @inheritParams sample_path_in_interval
#' @param init_time
#' @param final_time
#' @param init_state
#' @param final_state
#' @param irm_key
#'   
#' @return updated configuration matrix with a valid path (no update if constant
#'   path)
#'   
#' @export

sample_forward <- function(epimodel, subject, subj_ID, init_time, final_time, init_state, final_state, irm_key) {
          
          .valid_path         <- FALSE
          
          while(.valid_path == FALSE) {
                    
                    .keep_going         <- TRUE
                    .t                  <- init_time
                    .cur_state          <- init_state
                    .ind_cur            <- epimodel$.subj_row_ind
                    
                    while(.keep_going == TRUE) {
                              
                              if(.cur_state %in% epimodel$absorbing_states) {
                                        
                                        # stop forward sampling
                                        .keep_going = FALSE
                                        
                                        # determine if the path is valid. 
                                        if(.cur_state == final_state) {
                                                  .valid_path <- TRUE
                                                  epimodel$.subj_row_ind <- .ind_cur
                                                  
                                        } 
                                        
                                        break
                              } 
                              
                              # sample the next transition time
                              .t <- .t + rexp(1, rate = abs(epimodel$.irm[[irm_key]][.cur_state, .cur_state]))
                              
                              # if the next transition time is after the right endpoint, 
                              # stop sampling forward and determine if the path is valid
                              if(.t > final_time) {
                                        
                                        # stop forward sampling
                                        .keep_going = FALSE
                                        
                                        # determine if the path is valid. If valid, keep 
                                        # the extra rows and set epimodel$.subj_row_ind to 
                                        # .ind_cur. otherwise, .ind_cur and .cur_state will be 
                                        # reset in the next attempt
                                        if(.cur_state == final_state) {
                                                  .valid_path <- TRUE
                                                  epimodel$.subj_row_ind <- .ind_cur
                                        } 
                                        
                              # if .t <= final_time, sample the next state and 
                              # create an updated row in the config matrix
                              } else {
                                        # sample the next state
                                        .next_state <- sample.int(epimodel$num_states, 1, prob = pmax(epimodel$.irm[[irm_key]][.cur_state,], 0))
                                        
                                        # identify the event
                                        .event <- which((epimodel$flow[, .cur_state] == -1) & (epimodel$flow[, .next_state] == 1))
                                        
                                        # update the configuration matrix
                                        epimodel$config_mat[.ind_cur, c("time", "ID", "Event", subj_ID)] <- c(.t, subject, .event, .next_state)
                                        
                                        # update the sampling object
                                        .ind_cur            <- .ind_cur + 1
                                        .cur_state          <- .next_state 
                              }
                    }
          }
}