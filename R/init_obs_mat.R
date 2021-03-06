#' Initializes a bookkeeping object for the observed and augmented data.
#'
#' @param epimodel bookkeeping object.
#'
#' @return Bookkeeping object for observed and true counts at observation times.
#' @export
#'
init_obs_mat <- function(epimodel){

          # create object
          obs_mat <- matrix(0, nrow = length(epimodel$obstimes), ncol = (1 + 2*length(epimodel$meas_vars)))

          # insert obstimes
          obs_mat[,1] <- epimodel$obstimes

          # set column names
          colnames(obs_mat) <- c("time", paste(rep(epimodel$meas_vars, each = 2), c("_observed","_augmented"), sep = ""))

          # if a config_mat is provided, compute the true counts for the latent process at observation times
          if(!is.null(epimodel$pop_mat)) {
                    obs_mat[,paste0(epimodel$meas_vars,"_augmented")] <- 
                              epimodel$pop_mat[findInterval(epimodel$obstimes, epimodel$pop_mat[,"time"]),
                                               epimodel$meas_vars]
          }

          # if data was provided, insert it into the appropriate observed columns
          if(!is.null(epimodel$dat)) {
                    obs_mat[1:nrow(epimodel$dat), paste(epimodel$meas_vars, "_observed", sep="")] <- epimodel$dat[, epimodel$meas_vars]
          }

          return(obs_mat)
}