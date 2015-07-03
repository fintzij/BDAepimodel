#' Initializes a matrix of subject level trajectories.
#' 
#' @inheritParams simulate_epimodel
#'   
#' @return matrix of subject level trajectories, where each subject has two
#'   rows, one for time zero and one for tmax.
#' @export
#' 
#' @examples init_state <- c(S = 45, I = 5, R = 0)
#' tmax <- 5
#' init_subj_mat(init_state, tmax)
init_subj_mat <- function(init_state, tmax){
         
          # calculate popsize
          popsize <- sum(init_state)
          
          # set up matrix
          subj_mat <- matrix(0, nrow = 2*popsize, ncol = 2 + length(init_state))
          colnames(subj_mat) <- c("ID", "time", names(init_state))
          
          # individual IDs and times
          subj_mat[, 1] <- rep(1:popsize, each = 2)
          subj_mat[, 2] <- rep(c(0,tmax), popsize)
          
          # set initial and final statuses according to init_state
          IDs <- 1:popsize
          for(t in seq_along(init_state))  {
                    selected <- IDs[1:init_state[t]] 
                    subj_mat[subj_mat[ ,"ID"] %in% selected, names(init_state)[t]] <- 1
                    IDs <- IDs[-c(1:init_state[t])]
          }
}