#' Sample the infection status for the specified subject at observation times.
#'
#' @inheritParams draw_trajec
#'
#' @return updated configuration matrix
#' 
#' @export

sample_at_obs_times <- function(epimodel, subject, subj_ID) {
          
          epimodel$config_mat[epimodel$.ind_final_config, subj_ID] <- sample.int(n = epimodel$num_states, size = 1, prob = colSums(epimodel$.fb_mats[[epimodel$nobs - 1]]))
          
          for(s in (epimodel$nobs - 1) : 1) {
                    epimodel$config_mat[epimodel$.obs_time_inds[s], subj_ID] <- sample.int(n = epimodel$num_states, size = 1, prob = epimodel$.fb_mats[[s]][, epimodel$config_mat[epimodel$.obs_time_inds[s + 1], subj_ID]])
                    
          }
          
          # if the model is absorbing and progressive, fix the state in constant
          # inter-observation intervals and set the state after absorbtion
          if(epimodel$progressive & epimodel$absorbing_states) {
                    
                    # set the path to be constant in inter-observation intervals
                    # with identical endpoints
                    for(t in 1:(epimodel$nobs - 1)) {
                              
                              # check if the endpoints of the interval are equal
                              if(epimodel$config_mat[epimodel$.obs_time_inds[t], subj_ID] == epimodel$config_mat[epimodel$.obs_time_inds[t + 1], subj_ID]) {
                                        
                                        # if so, set the intermediate states to be equal
                                        epimodel$config_mat[epimodel$.obs_time_inds[t] : epimodel$.obs_time_inds[t+1], subj_ID] <- epimodel$config_mat[epimodel$.obs_time_inds[t], subj_ID]
                                        
                              }
                              
                    }
          
          # if the model merely has an absorbing state but is not progressive,
          # fix the state after absorbtion
          } else if(epimodel$absorbing_states & !epimodel$progressive) {
                    
                    # identify the first observation time when the subject is in
                    # the absorbing state, possibly never
                    .absorb_ind <- which(epimodel$config_mat[epimodel$.obs_time_inds, subj_ID] %in% epimodel$absorbing_states)[1]
                    
                    # if the subject was observed to be in an absorbing state,
                    # set the states thereafter
                    if(length(.absorb_ind) != 0) {
                              
                              epimodel$config_mat[epimodel$.obs_time_inds[.absorb_ind] : epimodel$.ind_final_config, subj_ID] <- epimodel$config_mat[epimodel$.obs_time_inds[.absorb_ind], subj_ID]
                              
                    }
                    
          }
}