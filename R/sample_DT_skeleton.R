#' Sample the state at event times from the discete time skeleton.
#'
#' This function only gets called when there is at least one transition event in
#' an inter-observation interval
#'
#' @inheritParams simulate_epimodel
#' @param path column from the config_mat corresponding to a subject path
#' @param subject subject ID, i.e. column index in config_mat
#' @param init_ind,final_ind index in the configuration matrix of the left/right
#'   endpoints of the observation interval
#'
#' @return updated configuration matrix
#'
#' @export
#'
sample_DT_skeleton <- function(path, epimodel, init_ind, final_ind) {
        
        # set the initial and final states
        .init_state         <- path[init_ind]
        .final_state        <- path[final_ind]
        
        # set up subseq vector
        .subseq_vec         <- c(.init_state, rep(0, final_ind - init_ind - 1), .final_state)
        .cur_state          <- .subseq_vec[1]
        
        .ind <- 1
        
        # sample the state at event times
        for(s in (init_ind + 1):(final_ind - 1)) {
                
                .ind <- .ind + 1
                
                .state_probs <- epimodel$tpms[.cur_state, , s-1] * epimodel$tpm_products[, .final_state, s] / epimodel$tpm_products[.cur_state, .final_state, s-1]
                
                # replace NaN values arising from 0/0 divisions
                .state_probs <- replace(.state_probs, is.nan(.state_probs), 0)
                
                .subseq_vec[.ind] <- .cur_state <- sample.int(n = epimodel$num_states, size = 1, replace = FALSE, prob = .state_probs)
        }
        
        path[init_ind:final_ind] <- .subseq_vec
        
        return(path)
}
