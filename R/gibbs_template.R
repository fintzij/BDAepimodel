#' Template for Gibbs kernel
#' 
#' Template for writing a MH transition kernel, example based on SIR model with
#' Gibbs updates for the sampling parameter.
#' 
#' @param epimodel
#'   
#' @return Template for Gibbs updates - example given for SIR model
#' @export
#' 
gibbs_template <- function(epimodel) {
          
          ############################################################################
          ### TEMPLATE FOR GIBBS UPDATES TO MEASUREMENT PARAMETER IN SIR MODEL.    ###
          ### COPY AND PASTE THIS TEMPLATE INTO A NEW SCRIPT AND EDIT ACCORDINGLY. ###
          ### DO NOT EDIT IN PLACE!!!!                                             ###
          ############################################################################
          
          # Get current parameters
          proposal          <- epimodel$params
          
          # draw new parameters - edit this line appropriately
          proposal["rho"]   <- rbeta(1, shape1 = 1 + sum(epimodel$obs_mat[,"I_observed"]), shape2 = 1 + sum(epimodel$obs_mat[,"I_augmented"] - epimodel$obs_mat[,"I_observed"]))
          
          # Calculate the measurement process likelihood. 
          obs_likelihood_new <- calc_obs_likelihood(epimodel, params = proposal, log = TRUE)
          
          # uncomment the following lines if resampling a rate parameter
          # epimodel <- build_new_irms(epimodel, proposal)
          # pop_likelihood_new  <- calc_pop_likelihood(epimodel, log = TRUE) 

          # update parameters and likelihoods - DO NOT DELETE THIS LINE 
          # if sampling a parameter that affects the calculation of the 
          #  population level likelihood, include that as well
          # epimodel <- update_params(epimodel, params = proposal, pop_likelihood = pop_likelihood_new, obs_likelihood = obs_likelihood_new)
          
          epimodel <- update_params(epimodel, params = proposal, obs_likelihood = obs_likelihood_new)
          
          return(epimodel)
}