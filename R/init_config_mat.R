#' initializes the bookkeeping object for population level counts
#'
#' @inheritParams init_epimodel
#' @param epimodel an epimodel list
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
init_config_mat <- function(epimodel, init_state, t0, tmax){

          # initialize matrix object
          pop_mat <- matrix(rep(c(t0, 0, 0, init_state), 2), nrow = 2, byrow = TRUE)
          config_mat <- matrix(as.integer(rep(rep(1:length(init_state), init_state), 2)), nrow = 2, byrow = TRUE)

          # set column names
          colnames(pop_mat) <- c("time", "ID", "Event", names(init_state))
          colnames(config_mat) <- paste0(".X", 1:sum(init_state))

          # set tmax
          pop_mat[2, "time"] <- tmax

          # pass to epimodel
          epimodel$pop_mat <- pop_mat
          epimodel$config_mat <- config_mat

          return(epimodel)
}