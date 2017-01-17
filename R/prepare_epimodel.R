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
          # observation times and instatiates epimodel$ind_final_config
          epimodel <- expand_pop_mat(epimodel)
          
          # initialize vector for storing the current subject level path
          epimodel$subj_path <- rep(0L, nrow(epimodel$pop_mat)) 
          
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
          epimodel <- build_new_irms(epimodel, epimodel$params)
          
          # update the eigen decompositions
          if(is.null(epimodel$sim_settings$analytic_eigen)) {
                    
                    buildEigenArray(real_eigenvals = epimodel$real_eigen_values,
                                    imag_eigenvals = epimodel$imag_eigen_values,
                                    eigenvecs      = epimodel$eigen_vectors, 
                                    inversevecs    = epimodel$inv_eigen_vectors, 
                                    irm_array      = epimodel$irm, 
                                    n_real_eigs    = epimodel$n_real_eigs)  
                    
          } else if(epimodel$sim_settings$analytic_eigen == "SIR") {
                    
                    buildEigenArray_SIR(real_eigenvals = epimodel$real_eigen_values,
                                        imag_eigenvals = epimodel$imag_eigen_values,
                                        eigenvecs      = epimodel$eigen_vectors, 
                                        inversevecs    = epimodel$inv_eigen_vectors, 
                                        irm_array      = epimodel$irm, 
                                        n_real_eigs    = epimodel$n_real_eigs,
                                        initial_calc   = TRUE)
                    
          } else if(epimodel$sim_settings$analytic_eigen == "SEIR") {
                    
                    buildEigenArray_SEIR(real_eigenvals = epimodel$real_eigen_values,
                                         imag_eigenvals = epimodel$imag_eigen_values,
                                         eigenvecs      = epimodel$eigen_vectors, 
                                         inversevecs    = epimodel$inv_eigen_vectors, 
                                         irm_array      = epimodel$irm, 
                                         n_real_eigs    = epimodel$n_real_eigs,
                                         initial_calc   = TRUE)
                    
          } 
          
          # get the irm keys
          epimodel$keys <- retrieveKeys(1:epimodel$ind_final_config, epimodel$irm_key_lookup, epimodel$pop_mat, epimodel$index_state_num)
          
          # compute the population level likelihood, and the measurement
          # process likelihood
          epimodel$likelihoods$pop_likelihood_cur <-
                    populationLikelihood(
                              pop_mat             = epimodel$pop_mat,
                              irm                 = epimodel$irm,
                              initdist            = epimodel$initdist,
                              initdist_param_inds = epimodel$initdist_param_inds,
                              flow_inds           = epimodel$flow_inds,
                              keys                = epimodel$keys,
                              inds                = c(1,
                                                      c(1:epimodel$ind_final_config)[-epimodel$obs_time_inds],
                                                      epimodel$ind_final_config),
                              loglik = TRUE
                    )
          
          epimodel$likelihoods$obs_likelihood <- calc_obs_likelihood(epimodel, log = TRUE)

          return(epimodel)
}