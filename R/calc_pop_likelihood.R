#' Calculate likelihood for the population level trajectory
#'
#' @inheritParams simulate_epimodel
#' @param params optional vector of parameters if values other than the current
#'   values should be used (e.g. in a transition kernel)
#' @param log TRUE/FALSE for whether to return the log-likelihood or likelihood.
#'
#' @return vector with log-likelihood of the population level trajectory and
#'   measurement process log-likelihood.
#'
#' @export
#'

calc_pop_likelihood <- function(epimodel, log = TRUE) {

         pop_likelihood <- populationLikelihood(pop_mat = epimodel$pop_mat, irm_array = epimodel$irm, initdist = epimodel$initdist, initdist_param_inds = epimodel$initdist_param_inds, flow_inds = epimodel$flow_inds, keys = epimodel$keys, inds = c(1, c(1:epimodel$ind_final_config)[-epimodel$obs_time_inds], epimodel$ind_final_config), loglik = log)

          return(pop_likelihood)
}