#' Calculate the likelihood for a subject-level trajectory based on a
#' time-inhomogeneous CTMC.
#' 
#' @inheritParams calc_pop_likelihood
#' @param subj_ID string for subject ID
#'   
#' @return likelihood or log-likelihood

calc_subj_likelihood <- function(epimodel, subject, subj_ID, log = TRUE) {
          
          # get the appropriate indices to exclude the intermediate observation time rows
          .left_endpoints     <- c(1, which(epimodel$config_mat[,"ID"] != 0))
          .right_endpoints    <- c(.left_endpoints[-1], epimodel$.ind_final_config)
          
          # get the subject_path
          .subj_path <- epimodel$config_mat[c(1, .right_endpoints), subj_ID]
          
          # compute the time diffs
          .time_diffs <- epimodel$config_mat[.right_endpoints, "time"] - epimodel$config_mat[.left_endpoints, "time"]
          
          # get irm keys
          .irm_keys <- generate_keys(epimodel, inds = .left_endpoints)
          
          # get hazards
          .hazards            <- rep(0, length(.left_endpoints))
          .event_rates        <- rep(0, length(.left_endpoints))
          
          for(r in 1:length(.hazards)) {
                    
                    # get the hazard from the appropriate rate matrix
                    .hazards[r] <- epimodel$.irm[[.irm_keys[r]]][.subj_path[r], .subj_path[r]]
                    
                    # check if the event corresponded to the subject. if so, get the event rate
                    if(epimodel$config_mat[.right_endpoints[r], "ID"] == subject) {
                              
                              .event_rates[r] <- epimodel$.irm[[.irm_keys[r]]][.subj_path[r], .subj_path[r + 1]]
                                        
                    }
          }
          
          
          subj_likelihood <- epimodel$d_initdist(state = .subj_path[1], params = epimodel$params, log = TRUE) + 
                    sum(ifelse(.event_rates == 0, .hazards * .time_diffs, log(.event_rates) + .hazards * .time_diffs))
          
          if(log == FALSE) {
                    subj_likelihood <- exp(subj_likelihood)
          }
          
          return(subj_likelihood)
}