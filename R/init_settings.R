#' Initialize a bookkeeping list of MCMC settings.
#' 
#' @param niter number of iterations/parameter updates for which to run the 
#'   sampler.
#' @param params_every thin parameter updates by saving every k-th draw.
#' @param configs_every thin configurations by saving the current configuration 
#'   after every k-th parameter update.
#' @param kernel transition kernel for parameters, provided as a list of 
#'   transition functions. Each function should take as an argument only the 
#'   epimodel list. The kernels should be self-contained, in particular making 
#'   any necessary transformations of parameters and calculating prior 
#'   probabilities as needed. Each kernel function should include a call to 
#'   \code{update_params} at the end
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
#' @param seed optional seed value. A random seed is generated and saved if none
#'  is supplied.
#'   
#' @return bookkeeping list
#' @export
#' 
#' 
init_settings <- function(niter, save_params_every = 1, save_configs_every = niter %/% 100, kernel, cov_mtx = NULL, configs_to_redraw, to_estimation_scale = NULL, from_estimation_scale = NULL, seed = NULL) {
          
          # check that if the one of the transformation arguments was provided,
          # then its inverse should be as well.
          if(is.null(to_estimation_scale) & !is.null(from_estimation_scale)){
                    stop(sQuote("to_estimation_scale"),"must be provided if from_estimation_scale is specified")    
          }
          
          if(is.null(from_estimation_scale) & !is.null(to_estimation_scale)){
                    stop(sQuote("from_estimation_scale"),"must be provided if to_estimation_scale is specified")
          }
          
          if(is.null(seed)) {
                    seed <- round(runif(1) * 1e8)
          }
          
          # initialize simulation settings list
          sim_settings <- list(niter = niter,
                               save_params_every = save_params_every,
                               save_configs_every = save_configs_every,
                               kernel = kernel,
                               cov_mtx = cov_mtx,
                               configs_to_redraw = configs_to_redraw,
                               to_estimation_scale = to_estimation_scale,
                               from_estimation_scale = from_estimation_scale,
                               seed = seed)
          
          return(sim_settings)
}