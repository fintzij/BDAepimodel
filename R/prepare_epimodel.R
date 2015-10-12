#' Prepares an epimodel environment. 
#'
#' @param epimodel initialized epimodel list
#'
#' @return epimodel environment
#' @export
#'
prepare_epimodel <- function(epimodel) {
          
          # convert epimodel list to environment to enable multiple assignment
          .epimodel <- list2env(epimodel, parent = emptyenv(), hash = TRUE)
          
          # initialize the list of log-likelihoods
          .epimodel$likelihoods <- list(pop_likelihood_cur = NULL,
                                        pop_likelihood_new = NULL,
                                        subj_likelihood_cur = NULL,
                                        subj_likelihood_new = NULL,
                                        obs_likelihood  = NULL)
          
          
          # identify whether there is structure in the flow matrix that can be
          # leveraged when resampling subject-level trajectories
          detect_structure(.epimodel)
          
          # add buffer to the configuration matrix, also adds configurations at
          # observation times and instatiates .epimodel$.ind_final_config
          .epimodel$config_mat <- expand_config_mat(.epimodel)
          
          # Initialize two lists for the transition probability matrices, one 
          # for the matrices and one for products. Also initialize a bookkeeping
          # vector indicating which tpms and tpm product matrices need to be 
          # recomputed. All tpms and products need to be recomputed following 
          # parameter updates, but only the tpms and products for the intervals 
          # affected by a change in one of the compartments indicated by 
          # index_states need to be recomputed when a new trajectory is redrawn.
          .epimodel$.tpms               <- vector("list", length = nrow(.epimodel$config_mat))
          .epimodel$.tpm_products       <- vector("list", length = nrow(.epimodel$config_mat))
          
          # initialize the matrix of emission probabilities
          .epimodel$nobs                <- length(.epimodel$obstimes)
          .epimodel$.emission_mat       <- matrix(0, nrow = .epimodel$num_states, ncol = .epimodel$nobs, dimnames = list(.epimodel$states, .epimodel$obstimes))
          
          # initialize the forward-backward matrices
          .epimodel$.fb_mats             <- vector("list", length = .epimodel$nobs - 1)
          
          # initialize the vector of subject-level initial state probabilities
          .epimodel$.initdist <- build_initdist(.epimodel)
          
          # get indices for the subject configuration portion of the configuration matrix
          .epimodel$.config_inds <- which(grepl(".X", colnames(.epimodel$config_mat)))
          
          # compute the population level likelihood, and the measurement
          # process likelihood
          .epimodel$likelihoods$pop_likelihood_cur <- calc_pop_likelihood(epimodel = .epimodel, log = TRUE)
          .epimodel$likelihoods$obs_likelihood <- calc_obs_likelihood(.epimodel, log = TRUE)
          
          
          return(.epimodel)
}