#' Prepares an epimodel environment.
#'
#' @param epimodel initialized epimodel list
#'
#' @return epimodel environment
#' @export
#'
prepare_epimodel <- function(epimodel) {

           # initialize the list of log-likelihoods
          epimodel$likelihoods <- list(pop_likelihood_cur = NULL,
                                        pop_likelihood_new = NULL,
                                        subj_likelihood_cur = NULL,
                                        subj_likelihood_new = NULL,
                                        obs_likelihood  = NULL)


          # identify whether there is structure in the flow matrix that can be
          # leveraged when resampling subject-level trajectories
          epimodel <- detect_structure(epimodel)

          # add buffer to the configuration matrix, also adds configurations at
          # observation times and instatiates epimodel$.ind_final_config
          epimodel <- expand_config_mat(epimodel)

          # Initialize two arrays for the transition probability matrices, one
          # for the matrices and one for products. Also initialize a bookkeeping
          # vector indicating which tpms and tpm product matrices need to be
          # recomputed. All tpms and products need to be recomputed following
          # parameter updates, but only the tpms and products for the intervals
          # affected by a change in one of the compartments indicated by
          # index_states need to be recomputed when a new trajectory is redrawn.
          epimodel$tpms               <- array(0, dim = c(epimodel$num_states, epimodel$num_states, nrow(epimodel$pop_mat)))
          epimodel$tpm_products       <- array(0, dim = c(epimodel$num_states, epimodel$num_states, nrow(epimodel$pop_mat)))

          # initialize the matrix of emission probabilities
          epimodel$nobs               <- length(epimodel$obstimes)
          epimodel$emission_mat       <- matrix(0, nrow = epimodel$num_states, ncol = epimodel$nobs, dimnames = list(epimodel$states, epimodel$obstimes))

          # initialize the forward-backward matrices
          epimodel$fb_mats             <- array(0, dim = c(epimodel$num_states, epimodel$num_states, epimodel$nobs - 1))

          # initialize the vector of subject-level initial state probabilities
          epimodel$initdist <- build_initdist(epimodel)

          # initialize objects for storing rate matrices and eigen decompositions
          epimodel <- init_irm(epimodel)
          
          # compute the rates
          rates <- do.call(cbind, lapply(epimodel$rates, do.call, list(state = epimodel$irm_key_lookup, params = epimodel$params)))
          
          # update the rate matrices
          buildRateArray(irm_array = epimodel$irm, rates = rates, flow_inds = epimodel$flow_inds)
          
          epimodel$keys <- retrieveKeys(1:epimodel$ind_final_config, epimodel$irm_key_lookup, epimodel$pop_mat, epimodel$index_state_num)
          
          # compute the population level likelihood, and the measurement
          # process likelihood
          epimodel$likelihoods$pop_likelihood_cur <- calc_pop_likelihood(epimodel = epimodel, log = TRUE)
          epimodel$likelihoods$obs_likelihood <- calc_obs_likelihood(epimodel, log = TRUE)

          return(epimodel)
}