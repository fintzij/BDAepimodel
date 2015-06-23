#' initializes the bookkeeping object for population level counts
#' 
#' @inheritParams simulate_epimodel
#'   
#' @return initialized matrix with columns to store event times and counts of 
#'   individuals in each compartment. The first and last rows of the matrix 
#'   always correspond to time 0 and tmax. The state at tmax is initialized to
#'   init_state.
#' @export
#' 
#' @examples init_state <- c(S = 45, I = 5, R = 0)
#' tmax <- 5
#' init_pop_traj(init_state, tmax)
#' 
init_pop_mat <- function(init_state, tmax){
          
          # initialize matrix object
          pop_mat <- matrix(rep(c(0, init_state), 2), nrow = 2, ncol = 1 + length(init_state), byrow = TRUE)
          pop_mat[2, 1] <- tmax
          
          # set column names
          colnames(pop_mat) <- c("time", names(init_state))
          
          return(pop_mat)
}