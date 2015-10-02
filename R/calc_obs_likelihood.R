#' Calculate likelihood for the population level trajectory and the measurement
#' process.
#' 
#' @inheritParams simulate_epimodel
#' @inheritParams calc_pop_likelihood
#'   
#' @return measurement process likelihood or log-likelihood.
#' 
#' @export
#' 
calc_obs_likelihood <- function(epimodel, log = FALSE) {
          
          obs_likelihood <- sum(epimodel$d_meas_process(state = epimodel$obs_mat, meas_vars = epimodel$meas_vars, params = epimodel$params, log = TRUE))
                    
          if(log == FALSE) {
                    obs_likelihood <- exp(obs_likelihood)
          }

          return(obs_likelihood)
                    
}