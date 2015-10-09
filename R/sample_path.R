#' Wrapper for sampling the path within homogeneous intervals. 
#'
#' @inheritParams draw_trajec
#' @param subj_ID character string for subject ID code
#'
#' @return updated configuration matrix within the epimodel environment. 
#' 
#' @export

sample_path <- function(epimodel, subject, subj_ID) {
         
          # expand the configuration matrix if the buffer is less than 10% of
          # the size of the current configuration
          if((nrow(epimodel$config_mat) - epimodel$.ind_final_config) < ceiling(0.1 * epimodel$.ind_final_config)) {
                    expand_config_mat(epimodel, buffer_size = ceiling(0.1 * epimodel$.ind_final_config))
          }
          
          
          # if the model is progressive and has an absorbing state, only need to
          # sample times of state transition in intervals with a state change
          if(epimodel$progressive & epimodel$absorbing_states) {
                    
                    for(s in which(diff(epimodel$config_mat[1:epimodel$.ind_final_config, subj_ID], lag = 1) != 0)) {
                              
                              sample_path_in_interval(epimodel, subject, subj_ID, interval = s)
                              
                    }
                    
          
          # if the model has an absorbing state, but is not progressive, no need
          # to sample times of transition after absorbtion
          } else if(epimodel$absorbing_states & !epimodel$progressive) {
                    
                    for(s in setdiff(1:(epimodel$.ind_final_config - 1), which(epimodel$config_mat[, subj_ID] %in% epimodel$absorbing_states)[1] : (epimodel$.ind_final_config - 1))) {
                              
                              sample_path_in_interval(epimodel, subject, subj_ID, interval = s)
                              
                    }
                    
                    
          # if the model does not have an absorbing state, sample paths in all
          # inter-event intervals
          } else {
                    
                    for(s in 1:(epimodel$.ind_final_config - 1)) {
                              
                              sample_path_in_interval(epimodel, subject, subj_ID, interval = s)
                              
                    }
                    
          }
          
}