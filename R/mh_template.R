#' Template for Metropolis-Hastings kernel
#' 
#' Template for writing a MH transition kernel, example based on SIR model with
#' MH proposals for infectivity and recovery parameters on log transformed scale.
#' 
#' @param epimodel
#'   
#' @return template for Metropolis-Hastings updates - example given for SIR model
#' @export
#' 
mh_template <- function(epimodel) {
          
          ############################################################################
          ### TEMPLATE FOR METROPOLIS-HASTINGS PROPOSALS IN SIR MODEL.             ###
          ### COPY AND PASTE THIS TEMPLATE INTO A NEW SCRIPT AND EDIT ACCORDINGLY. ###
          ### DO NOT EDIT IN PLACE!!!!                                             ###
          ############################################################################
          
          # get existing parameter values 
          proposal <- epimodel$params
          
          # get prior likelihoods - edit this line appropriately
          prior_prob_cur <- c(dnorm(epimodel$params["beta"], mean = log(1), sd = 1.8, log = TRUE),
                              dnorm(epimodel$params["mu"], mean = log(1), sd = 3, log = TRUE)) 
          
          # set log likelihood function - DO NOT DELETE THIS LINE
          log_likelihood_cur <- epimodel$likelihoods$pop_likelihood_cur + epimodel$likelihoods$obs_likelihood
          
          # transform to estimation scale -- comment out if no transformation needed
          proposal <- epimodel$sim_settings$to_estimation_scale(proposal, epimodel)
          
          # make proposal - edit this line appropriately
          proposal[c("beta", "mu")] <- MASS::mvrnorm(n = 1, mu = proposal[c("beta", "mu")], Sigma = epimodel$sim_settings$cov_mtx)
          
          # get prior likelihood for proposed parameters - edit appropriately
          prior_prob_new <- c(dnorm(proposal[[1]], mean = log(1), sd = 1.8, log = TRUE),
                              dnorm(proposal[[2]], mean = log(1), sd = 3, log = TRUE))
          
          # transform back to the model scale - comment out if no transformation needed
          proposal <- epimodel$sim_settings$from_estimation_scale(proposal, epimodel)
          
          
          # update array of rate matrices - DO NOT DELETE THIS LINE IF ANY RATE PARAMETERS HAVE BEEN INCREMENTED
          epimodel <- build_new_irms(epimodel, proposal)
          
          # get population level complete data likelihood for the new parameters
          # if the population-level CTMC likelihood or the measurement process
          # is unchanged, comment out the corresponding line
          pop_likelihood_new  <- calc_pop_likelihood(epimodel, log = TRUE) 
          obs_likelihood_new  <- calc_obs_likelihood(epimodel = epimodel, params = proposal, log = TRUE)
          log_likelihood_new  <- pop_likelihood_new + obs_likelihood_new
          
          # compute the acceptance probability and accept or reject proposal 
          accept_prob <- log_likelihood_new - log_likelihood_cur + sum(prior_prob_new) - sum(prior_prob_cur)
          
          if(accept_prob >= 0 || accept_prob >= log(runif(1))) {
                    
                    # update parameters and likelihood objects - no need to 
                    # include unchanged observation or population likelihoods
                    epimodel <- update_params(epimodel, params = proposal, pop_likelihood = pop_likelihood_new, obs_likelihood = obs_likelihood_new)
                    
                    # update the eigen decompositions - comment out this line if
                    # the proposal parameters are unchanged
                    buildEigenArray(eigenvals = epimodel$eigen_values, eigenvecs = epimodel$eigen_vectors, inversevecs = epimodel$inv_eigen_vectors, irm_array = epimodel$irm)
                    
          } else {
                    # restore the previous array of rate matrices - comment out
                    # this line if the proposal leaves the rate parameters
                    # unchanged
                    epimodel$irm <- epimodel$.irm_cur
          }
          
          return(epimodel)
}