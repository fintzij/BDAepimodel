#' Initializes a bookkeeping object for the observed and augmented data. 
#'
#' @param epimodel bookkeeping object. 
#'
#' @return Bookkeeping object for observed and true counts at observation times. 
#' 
.init_obs_mat <- function(epimodel){
          
          # create object
          obs_mat <- matrix(0, nrow = length(epimodel$dat[,"obstimes"]), ncol = (1 + 2*length(epimodel$meas_vars)))
          
          # insert obstimes
          obs_mat[,1] <- epimodel$dat[,"obstimes"]
          
          # set column names
          colnames(obs_mat) <- c("time", paste(rep(epimodel$meas_vars, each = 2), c("_observed","_truth"), sep = ""))
          
          # if data was provided, insert it into the appropriate observed columns
          obs_mat[, paste(epimodel$meas_vars, "_observed", sep="")] <- epimodel$dat[, meas_vars]
          
          return(obs_mat)
}