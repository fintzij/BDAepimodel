#' Sample the infection status at times of state transition of other subjects
#'
#' @inheritParams draw_trajec
#'
#' @return updated configuration matrix of subject level trajectories

sample_at_event_times <- function(epimodel, subject) {
          
          # get subject ID
          .subj_ID <- paste0(".X", subject)
          
          # If the model is progressive and has an absorbing state, only need to
          # sample the discrete time skeleton where the endpoints of observation
          # intervals indicate a jump
          if(epimodel$absorbing_states & epimodel$progressive) {
                    
                    for(t in 1:(epimodel$nobs - 1)) {
                              
                              # draw the status at inter-event times conditional on the
                              # status at observation times only if a change is observed
                              if(epimodel$config_mat[epimodel$.obs_time_inds[t], .subj_ID] != epimodel$config_mat[epimodel$.obs_time_inds[t + 1], .subj_ID]) {
                                        
                                        epimodel$config_mat[epimodel$.obs_time_inds[t] : epimodel$.obs_time_inds[t + 1], .subj_ID] <- sample_dt_skeleton(epimodel, init_ind = epimodel$.obs_time_inds[t], final_ind= epimodel$.obs_time_inds[t+1])
                              }
                              
                    }
          
          #if the model is not progressive but has an absorbing state,
          #no need to deal with the path after absorbtion
          } else if(epimodel$absorbing_states & !epimodel$progressive) {
                    
                    for(t in 1:(epimodel$nobs - 1)) {
                              
                              # draw the status at inter-event times conditional
                              # on the status at observation times only if the
                              # subject is not in an absorbing state
                              if(!epimodel$config_mat[epimodel$.obs_time_inds[t], .subj_ID] %in% epimodel$absorbing_states) {
                                        
                                        epimodel$config_mat[epimodel$.obs_time_inds[t] : epimodel$.obs_time_inds[t + 1], .subj_ID] <- sample_dt_skeleton(epimodel, init_ind = epimodel$.obs_time_inds[t], final_ind= epimodel$.obs_time_inds[t+1])
                                        
                              } else break
                              
                    }
                    
                    # if the model is both not progressive and lacks an
                    # absorbing state, need to sample all inter-event intervals
          } else {
                    for(t in 1:(epimodel$nobs - 1)) {
                              
                              epimodel$config_mat[epimodel$.obs_time_inds[t] : epimodel$.obs_time_inds[t + 1], .subj_ID] <- sample_dt_skeleton(epimodel, init_ind = epimodel$.obs_time_inds[t], final_ind= epimodel$.obs_time_inds[t+1])
                              
                    }
                    
          }
          
}