#' Calculate likelihood for the population level trajectory and the measurement process. 
#'
#' @inheritParams simulate_epimodel
#'
#' @return vector with log-likelihood of the population level trajectory and measurement process log-likelihood. 
#' @export
#'
calc_log_likelihoods <- function(epimodel, pop = TRUE, obs = TRUE) {
          
          if(pop == TRUE) {
                    epimodel$pop_log_likelihood 
          }
}