#' Wrapper for sampling the path within homogeneous intervals. 
#'
#' @inheritParams draw_trajec
#' @param subject character string for subject ID code
#'
#' @return updated configuration matrix within the epimodel environment. 
#' 
#' @export

sample_path <- function(epimodel, subject) {

          # if the model is progressive and has an absorbing state, only need to
          # sample times of state transition in intervals with a state change
          if(epimodel$progressive & epimodel$absorbing_states) {
                    
                    intervals <- which(diff(epimodel$subj_path[1:epimodel$ind_final_config], lag = 1) != 0)
                    path      <- vector(mode = "list", length = length(intervals))
                    
                    for(s in seq_along(path)) {
                              
                              path[[s]] <- sample_path_in_interval(epimodel, subject, interval = intervals[s])
                              
                    }
                    
          
          # if the model has an absorbing state, but is not progressive, no need
          # to sample times of transition after absorbtion
          } else if(epimodel$absorbing_states & !epimodel$progressive) {
                    
                    intervals  <- setdiff(1:(epimodel$ind_final_config - 1), which(epimodel$subj_path[1:epimodel$ind_final_config] %in% epimodel$absorbing_states)[1] : (epimodel$ind_final_config - 1))
                    path <- vector(mode = "list", length = length(intervals))
                    
                    for(s in seq_along(path)) {
                              
                              path[[s]] <- sample_path_in_interval(epimodel, subject, interval = intervals[s])
                              
                    }
                    
                    
          # if the model does not have an absorbing state, sample paths in all
          # inter-event intervals
          } else {
                    path <- vector(mode = "list", length = length(intervals))
                    
                    for(s in 1:(epimodel$ind_final_config - 1)) {
                              
                              path[[s]] <- sample_path_in_interval(epimodel, subject, interval = s)
                              
                    }
                    
          }
          
          if(!is.null(path)) {
                    path <- do.call(rbind, path)
          }
          
          return(path)
}