#' Initializes a bookkeeping object for the data
#'
#' @inheritParams simulate_epimodel
#'
#' @return Bookkeeping object for observed and true counts at observation times. 
#' @export
#'
#' @examples obstimes <- seq(0, 5, by = 0.5)
#' meas_vars <- "I"
#' init_data_obj(obstimes, meas_vars)
init_data_mat <- function(obstimes, meas_vars){
          
          # create object
          obs_mat <- matrix(0, nrow = length(obstimes), ncol = (1 + 2*length(meas_vars)))
          
          # insert obstimes
          obs_mat[,1] <- obstimes
          
          # set column names
          colnames(obs_mat) <- c("time", paste(rep(meas_vars, each = 2), c("_observed","_truth"), sep = ""))
          
          return(obs_mat)
}