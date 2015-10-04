#' Initialize a bookkeeping list of MCMC settings.
#' 
#' @param niter number of iterations/parameter updates for which to run the 
#'   sampler.
#' @param burnin length of burn-in run in terms of number of parameter updates.
#' @param params_every thin parameter updates by saving every k-th draw.
#' @param configs_every thin configurations by saving the current configuration 
#'   after every k-th parameter update.
#' @param kernel transition kernel for parameters, provided as a list of 
#'   transition functions. Each function should use objects that are available 
#'   within an epimodel list or that are supplied as arguments to 
#'   init_simulation.
#' @param priors list of functions to evaluate prior distributions.
#' @param cov_mtx covariance matrix for parameters to be used in M-H updates.
#' @param configs_to_redraw number of subject level trajectories to sample 
#'   between parameter updates.
#' @param to_estimation_scale list of functions for transforming model 
#'   parameters to the scale on which new parameters should be proposed. List 
#'   element names should correspond to exactly to the parameter names given in 
#'   \code{params}.
#' @param from_estimation_scale list of functions for transforming model 
#'   parameters from the a scale on which new parameters are proposed to the 
#'   scale used to evaluate the process likelihood. List element names should 
#'   correspond to exactly to the parameter names given in \code{params}.
#' @param parallelize option to parallelize computation of transition 
#'   probability and rate matrices. Defaults to FALSE.
#'   
#' @return bookkeeping list
#' @export
#' 
#' @examples  
#' # Transition kernel for the SIR model - block MH updates for beta and mu, 
#' # Gibbs updates for rho and p0. cov_mtx is supplied as an additional argument
#' # to init_simulation. This function uses the mvrnorm function from the MASS 
#' # package.
#' 
#' # functions to go to and from the estimation scale
# to_estimation_scale <- list(function(params, popsize) {log(params["beta"] *
# popsize/params["mu"])}, function(params) {log(epimodel$params["mu"])})
# 
# from_estimation_scale <- list(function(params_est, popsize)
# {exp(params_est["mu"])/popsize * exp(params_est["beta"])}, 
# function(params_est) {exp(params_est["mu"])})
# 
# 
# beta_mu_kernel <- function(epimodel, cov_mtx, log_likelihood_cur, 
# prior_prob_cur, dprior, to_estimation_scale, from_estimation_scale) {
#        
#        # four vectors - existing parameters on natural and estimation scale, 
#        # and new parameters on natural and estimation scale. 
#        params_est <- epimodel$params
#        params_est[c("beta","mu")] <- to_estimation_scale(params_est[c("beta",
#                                                 "mu")], epimodel$popsize) 
#                                                 
#        proposal_est <- proposal <- params_est
#        
#        # make proposal
#        proposal_est[c("beta", "mu")] <- mvrnorm(n = 1, mu = params[c("beta", 
#        "mu")], Sigma = cov_mtx)
#        
#        proposal[c("beta","mu")]<-from_estimation_scale(proposal_est[c("beta",
#                                                 "mu")], epimodel$popsize)
#        
#        # update prior probability on estimation scale
#        prior_prob_new <- update_prior_probs(params_new = proposal_est, 
#                  params_cur = params_est, dprior = dprior, log = TRUE)
#        log_likelihood_new <- calc_log_likelihood(epimodel = epimodel, params 
#        = proposal)
#        
#        # compute the acceptance probability and accept or reject proposal
#        accept_prob <- log_likelihood_new - log_likelihood_cur + 
#        prior_prob_new - prior_prob_cur
#        
#        if(accept_prob >= 0 || accept_prob >= log(runif(1))) {
#                  params_new <- proposal
#        } else {
#                  params_new <- epimodel$params
#        }
#        
#        return(params_new)        
# }
# 
# rho_p0_kernel <- function(epimodel) {
#        
#        # Gibbs updates for rho and p0 - beta(1, 1) prior for each
#        params_new <- epimodel$params
#        params_new["rho"] <- rbeta(1, shape1 = 1 + sum(epimodel$obs_mat
#        [,"I_observed"]), shape2 = 1 + sum(epimodel$obs_mat[,"I_truth"] - 
#        epimodel$obs_mat[,"I_observed"]))
#        params_new["p0"] <- rbeta(1, shape1 = 1 + epimodel$config_mat[1, "I"],
#         shape2 = 1 + epimodel$popsize - epimodel$config_mat[1, "I"])
#        
#        return(params_new)        
# }
# 
# kernel <- list(beta_mu_kernel, rho_p0_kernel)
# 
# configs_to_redraw <- 10
#  
# # save new parameters every iteration
# # save every tenth configuration matrix
# # resample 10 subject-level trajectories in between parameter updates
# init_simulation(niter = 1000,
#                  burnin <- 200,
#                  params_every <- 1, 
#                  configs_every <- 10,
#                  kernel = kernel,
#                  configs_to_redraw = 10)
#' 
init_settings <- function(epimodel, niter, burnin = NULL, params_every = 1, configs_every = floor(niter / 100), kernel, cov_mtx = NULL, configs_to_redraw, to_estimation_scale = NULL, from_estimation_scale = NULL) {
          
          # check that the appropriate argumentsare provided as integers
          if(!all(c(niter, burnin, params_every, configs_every, configs_to_redraw) == floor(c(niter, burnin, params_every, configs_every, configs_to_redraw)))) {
                    which_not_int <- which(c(niter, burnin, params_every, configs_every, configs_to_redraw) != floor(c(niter, burnin, params_every, configs_every, configs_to_redraw)))
                    stop(paste(c(niter, burnin, params_every, configs_every, configs_to_redraw)[which_not_int], collapse = ", "), "must be integer valued.")
          }
          
          # check that if the one of the transformation arguments was provided,
          # then its inverse should be as well.
          if(is.null(to_estimation_scale) & !is.null(from_estimation_scale)){
                    stop(sQuote("to_estimation_scale"),"must be provided if from_estimation_scale is specified")    
          }
          
          if(is.null(from_estimation_scale) & !is.null(to_estimation_scale)){
                    stop(sQuote("from_estimation_scale"),"must be provided if to_estimation_scale is specified")
          }
          
          # initialize simulation settings list
          sim_settings <- list(niter = niter,
                               burnin = burnin, 
                               params_every = params_every,
                               configs_every = configs_every,
                               kernel = kernel,
                               cov_mtx = cov_mtx,
                               configs_to_redraw = configs_to_redraw,
                               to_estimation_scale = to_estimation_scale,
                               from_estimation_scale = from_estimation_scale)
          
          return(sim_settings)
}