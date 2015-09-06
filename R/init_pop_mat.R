#' initializes the bookkeeping object for population level counts
#' 
#' @inheritParams init_epimodel
#' @param t0 the first observation time
#' @param tmax the final observation time
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
init_config_mat <- function(init_state, t0, tmax){
          
          # initialize matrix object
          pop_mat <- matrix(rep(c(t0, 0, init_state, rep(1:length(init_state), init_state)), 2 * sum(init_state)), nrow = 2 * sum(init_state), byrow = TRUE)
          
          # set column names
          colnames(pop_mat) <- c("time", "ID", names(init_state), paste("X", 1:sum(init_state), sep = ""))
          # set tmax
          pop_mat[(sum(init_state) + 1):nrow(pop_mat), "time"] <- tmax
          
          # set IDs
          pop_mat[, "ID"] <- rep(1:sum(init_state), 2)
          
          return(pop_mat)
}