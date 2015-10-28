#' Initialize tuning for M-H kernels.
#' 
#' @param epimodel
#' @param tune_params character vector of parameters for which proposal is to be
#'   tuned.
#' @param cov_mtx covariance matrix for parameters to be used in M-H updates - 
#'   supercedes the covariance matrix that already exists in the epimodel object
#'   if it was specified in both.
#' @param cov_scale scale tuning parameter for the covariance matrix - i.e. c * 
#'   Sigma, defaults to (2.38^2 / length(tune_params)).
#' @param mintune minimum number of tuning loops, defaults to 2.
#' @param maxtune maximum number of tuning loops, defaults to 24.
#' @param ntu number of parameter updates per tuning loop, defaults to 500.
#' @param tunewt initial weight to put on covariance matrix from most recent 
#'   tuning loop if it is desired to tune the covariance matrix according to a 
#'   weighted average of the previous covariance matrix and the new one. 
#'   Defaults to 0 - i.e. scale tuning only. The tuning weight will decrease in 
#'   proportion to the total number of iterations run - e.g. if the initial
#'   weight is 1, the covariance matrix will receive weight 1. The next 500 will
#'   receive weight 0.5. the 500 after that will receive weight 500/1500.
#' @param target_accept vector of length 2 giving the acceptable range for the 
#'   acceptance rate for M-H parameter updates, defaults to c(0.15, 0.5). Tuning
#'   ends once the acceptance rate either falls within this range or the maximum
#'   number of tuning loops is exhausted.
#'   
#' @return epimodel object with initialized tuning settings
#' @export
#'   
#' @examples 
#' \dontrun{init_tuning(epimodel, 
#'                  tune_params = c("beta", "mu", "rho"),
#'                  cov_mtx = matrix(c(1, -0.865, -0.747, 
#'                                      -0.865, 1, 0.831,
#'                                      -0.747, 0.831, 1), nrow = 3),
#'                  cov_scale = 0.02, 
#'                  mintune = 2,
#'                  maxtune = 10, 
#'                  ntu = 500,
#'                  tunewt = 0,
#'                  target_accept = 0.15, 0.5)}

init_tuning <- function(epimodel, tune_params, cov_mtx = NULL, cov_scale = (2.38^2 / length(tune_params)), mintune = 2, maxtune = 24, ntu = 500, tunewt = 0, target_accept = c(0.15, 0.5)) {
          
          if(is.null(cov_mtx) && is.null(epimodel$sim_settings$cov_mtx)) {
                    stop(sQuote("cov_mtx"), "must be specified either within the tuning settings or within the simulation settings")
          }
          
          if(is.null(epimodel$sim_settings)) {
                    stop(sQuote("set the simulation settings prior to setting the tuning settings"))
          }
          
          # set the tuning pointer to TRUE
          epimodel$tune <- TRUE
          
          # set the covariance matrix if one is specified
          if(!is.null(cov_mtx)) {
                    epimodel$sim_settings$cov_mtx <- cov_scale^2 * cov_mtx
          }
          
          epimodel$tune_settings <- list(tune_params = tune_params,
                                         cov_scale = cov_scale, 
                                         mintune = mintune,
                                         maxtune = maxtune,
                                         ntu = ntu,
                                         tunewt = tunewt,
                                         target_accept = target_accept)
          
          return(epimodel)
          
}