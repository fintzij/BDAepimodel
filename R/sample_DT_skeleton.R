#' Sample the state at event times from the discete time skeleton.
#' 
#' This function only gets called when there is at least one transition event in
#' an inter-observation interval
#' 
#' @inheritParams simulate_epimodel
#' @param subj_ID string indicating the subject ID, e.g. ".X1"
#' @param init_ind,final_ind index in the configuration matrix of the left/right
#'   endpoints of the observation interval
#'   
#' @return updated configuration matrix
#'   
sample_DT_skeleton <- function(epimodel, subj_ID, init_ind, final_ind) {
          
          # set the initial and final states
          .init_state         <- epimodel$config_mat[init_ind, subj_ID]
          .final_state        <- epimodel$config_mat[final_ind, subj_ID]
          
          # sample the state at event times
          for(s in (init_ind + 1):(final_ind - 1)) {
                    
                    .state_probs <- epimodel$.tpms[[s-1]][epimodel$config_mat[s-1, subj_ID], ] * epimodel$.tpm_products[[s]][, .final_state] / epimodel$.tpm_products[[s-1]][epimodel$config_mat[s-1, subj_ID], .final_state]
                    
                    # replace NaN values arising from 0/0 divisions
                    .state_probs <- replace(.state_probs, is.nan(.state_probs), 0)
                    
                    epimodel$config_mat[s, subj_ID] <- sample.int(n = epimodel$num_states, size = 1, prob = .state_probs)
          }
          
} 