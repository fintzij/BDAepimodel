#' Instatiates a data matrix based on a configuration matrix. The code in this
# function is duplicated in \code{init_obs_mat}, which was tested explicitly. 
#'
#' @inheritParams simulate_epimodel 
#' 
#' @return data matrix for compartment counts at observation times.
#'
#' @export
#'
traj_to_obs <- function(epimodel = NULL) {
          
          # ensure that a configuration matrix is supplied
          if(is.null(epimodel$config_mat)){
                    stop(sQuote("epimodel$config_mat"), "must be supplied")
          }
          
          # instatiate dat matrix
          dat <- matrix(0, nrow = length(epimodel$obstimes), ncol = 1 + length(epimodel$meas_vars))
          colnames(dat) <- c(epimodel$time_var, epimodel$meas_vars)
          
          dat[, epimodel$time_var] <- epimodel$obstimes
          
          # find the compartment count just before each observation time
          for(k in 1:length(obstimes)) {
                    dat[k, meas_vars] <- epimodel$config_mat[, meas_vars][sum(epimodel$config_mat[,"time"] <= dat[k, epimodel$time_var])]
          }
          
          return(dat)
}