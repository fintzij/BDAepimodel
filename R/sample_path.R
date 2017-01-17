#' Wrapper for sampling the path within homogeneous intervals. 
#'
#' @inheritParams draw_trajec
#' @param subject character string for subject ID code
#'
#' @return updated configuration matrix within the epimodel environment. 
#' 
#' @export

sample_path <- function(epimodel) {

          if(epimodel$absorbing_states) {
                    
                    # instatiate the vector intervals where the subject path needs 
                    # to be sampled and the matrix for storing the jumps
                    if(epimodel$progressive) {
                              
                              intervals    <- which(diff(epimodel$subj_path[1:epimodel$ind_final_config], lag = 1) != 0)
                              subject_path <- matrix(0.0, nrow = epimodel$num_states, ncol = 3)
                              
                    } else  {
                              intervals <- setdiff(1:(epimodel$ind_final_config - 1),
                                                  which(epimodel$subj_path[1:epimodel$ind_final_config] %in% 
                                                                  epimodel$absorbing_states)[1]:(epimodel$ind_final_config - 1))
                              subject_path <- matrix(0.0, nrow = epimodel$num_states, ncol = 3)
                    } 
          } else {
                    intervals    <- seq(1, epimodel$ind_final_config - 1)
                    subject_path <- matrix(0.0, nrow = 3 * epimodel$num_states, ncol = 3)
          }
                    
          # initial row index
          row_ind   <- 0
          n_jumps   <- 0
          path_rows <- nrow(subject_path)

          # sample the path in each interval 
          for(s in intervals) {
                    
                    init_state  <- epimodel$subj_path[s]
                    final_state <- epimodel$subj_path[s + 1]
                    
                    init_time   <- epimodel$pop_mat[s, "time"]
                    final_time  <- epimodel$pop_mat[s + 1, "time"]
                    
                    irm_key     <- epimodel$keys[s]
                    
                    # sample the path within the interval conditional on the state at the endpoints
                    if(epimodel$sim_settings$ecctmc_method == "mr") {
                              path <- ECctmc::sample_path_mr(init_state,
                                                             final_state,
                                                             init_time,
                                                             final_time,
                                                             epimodel$irm[, , irm_key])
                    } else {
                              path <- ECctmc::sample_path_unif3(init_state,
                                                                final_state,
                                                                init_time,
                                                                final_time,
                                                                epimodel$irm[, , irm_key],
                                                                epimodel$tpms[, , s])
                    }
                    
                    # compute the number of jumps
                    n_jumps <- nrow(path) - 2
                    
                    # insert jumps into the subject path matrix
                    if(n_jumps) {
                              
                              # add rows if necessary
                              if(row_ind + n_jumps > nrow(subject_path)) {
                                        subject_path <- rbind(subject_path, matrix(0.0, max(n_jumps, 2*epimodel$num_states), 3))
                              }
                              
                              # get indices for the copy
                              to_inds   <- row_ind + seq_len(n_jumps) # indices in the subject_path matrix
                              from_inds <- seq_len(n_jumps) + 1       # indices in the path matrix
                              
                              # insert the transition times and states into the sample path
                              subject_path[to_inds, c(1,3)] <- path[from_inds,]
                              
                              # now identify the event codes
                              for(t in seq_len(n_jumps)) {
                                        init_state <- path[from_inds[t] - 1, 2]
                                        next_state <- path[from_inds[t], 2]
                                        subject_path[to_inds[t], 2] <- which((epimodel$flow[, init_state] == -1) & 
                                                                                       (epimodel$flow[, next_state] == 1))
                              }
                              
                              # increment the row index
                              row_ind <- row_ind + n_jumps
                    }
          }
          
          return(subject_path[seq_len(row_ind),,drop = FALSE])
}