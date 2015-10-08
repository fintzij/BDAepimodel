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
#' @export
#'   
sample_DT_skeleton <- function(epimodel, subj_ID, init_ind, final_ind) {
          
          # set the initial and final states
          .init_state         <- epimodel$config_mat[init_ind, subj_ID]
          .final_state        <- epimodel$config_mat[final_ind, subj_ID]
          
          # set up subseq vector
          .subseq_vec         <- c(.init_state, rep(0, final_ind - init_ind - 1), .final_state)
          .cur_state          <- .subseq_vec[1]
          
          .ind <- 1

          # sample the state at event times
          for(s in (init_ind + 1):(final_ind - 1)) {
                    
                    .ind <- .ind + 1
                    
                    .one_step_tpm  <- epimodel$.tpms[[s-1]]
                    .tpm_prods_den <- epimodel$.tpm_products[[s-1]]
                    .tpm_prods_num <- epimodel$.tpm_products[[s]]
                    
                    .state_probs <- .one_step_tpm[.cur_state, ] * .tpm_prods_num[, .final_state] / .tpm_prods_den[.cur_state, .final_state]
                    
                    # replace NaN values arising from 0/0 divisions
                    .state_probs <- replace(.state_probs, is.nan(.state_probs), 0)
                    
                    .subseq_vec[.ind] <- .cur_state <- .Internal(sample(n = epimodel$num_states, size = 1, replace = FALSE, prob = .state_probs))
          }
          
          epimodel$config_mat[init_ind:final_ind, subj_ID] <- .subseq_vec
          
} 