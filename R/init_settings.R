#' Initialize a bookkeeping list of MCMC settings.
#' 
#' @param epimodel epimodel object
#' @param niter number of iterations/parameter updates for which to run the 
#'   sampler.
#' @param configs_to_redraw number of subject level trajectories to sample 
#'   between parameter updates.
#' @param preferential_sampling vector with two elements, given in the following
#'   order. First the size of the group to be preferentially sampled. Second, 
#'   the probability that a subject is sampled from the preferential group. The 
#'   IDs of subjects to be preferentially sampled will be selected uniformly at 
#'   random from among subjects with non-constant paths when the collection of 
#'   subject paths is initialized. If the size of the preferentially sampled 
#'   group is larger than the number of subjects with non-constant paths, 
#'   additional subject IDs will be added by selecting subejct IDs uniformly at 
#'   random from among the other subjects.
#' @param init_popsize population size with which to initialize the epidemic.
#' @param compartment_dist distribution according to which to place the 
#'   remainder of the subjects. required if init_popsize is not equal to the 
#'   true population size. Supplied as a named vector of probabilities.
#' @param save_params_every thin parameter updates by saving every k-th draw, 
#'   defaults to 1.
#' @param save_configs_every thin configurations by saving the current 
#'   configuration after every k-th parameter update.
#' @param kernel transition kernel for parameters, provided as a list of 
#'   transition functions. Each function should take as an argument only the 
#'   epimodel list. The kernels should be self-contained, in particular making 
#'   any necessary transformations of parameters and calculating prior 
#'   probabilities as needed. Each kernel function should include a call to 
#'   \code{update_params} at the end.
#' @param post_init_params a vector of initial parameter values to be used after
#'   an acceptable initial configuration has been simulated.
#' @param cov_mtx covariance matrix to be used in the transition kernels.
#' @param to_estimation_scale list of functions for transforming model 
#'   parameters to the scale on which new parameters should be proposed. List 
#'   element names should correspond to exactly to the parameter names given in 
#'   \code{params}.
#' @param from_estimation_scale list of functions for transforming model 
#'   parameters from the a scale on which new parameters are proposed to the 
#'   scale used to evaluate the process likelihood. List element names should 
#'   correspond to exactly to the parameter names given in \code{params}.
#' @param analytic_eigen optional. If NULL, eigenvalues, eigenvectors, and 
#'   inverses of the eigenvector matrices of each rate transition rate matrix 
#'   are computed numerically. If one of "SIR", or "SEIR", the eigen 
#'   decompositions will be computed analytically. It is up to the user to 
#'   ensure that the rate matrices are structured in the appropriate form (see 
#'   vignettes for examples).
#' @param ecctmc_method Method for sampling the exact times of state transition.
#'   Defaults to "mr" for modified rejection sampling. May also be specified as
#'   "unif" for uniformization.
#' @param seed optional seed value. A random seed is generated and saved if none
#'   is supplied.
#' @param time_limit optional maximum time budget for the MCMC, checked every 
#'   time the parameters are saved. Should be provided as a diffitme object with
#'   units in hours.
#'   
#' @return bookkeeping list
#' @export
#' 
init_settings <- function(epimodel, niter, configs_to_redraw, config_replacement = FALSE, preferential_sampling = NULL, init_popsize = NULL, compartment_dist = NULL, save_params_every = 1, save_configs_every = niter %/% 100, kernel, post_init_params = NULL, cov_mtx = NULL, to_estimation_scale = NULL, from_estimation_scale = NULL, analytic_eigen = NULL, ecctmc_method = "mr", seed = NULL, time_limit = NULL) {

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

          if(!is.null(preferential_sampling) && preferential_sampling[[1]] > epimodel$popsize) {
                    stop("The size of the group to be preferentially sampled cannot be larger than the population size.")
          }

          if(!is.null(init_popsize)) {
                    if(init_popsize > epimodel$pop) {
                              stop("The population size used to initialize the collection of paths must be no greater than the true population size.")
                    }

                    if(is.null(compartment_dist)) {
                              stop("The compartment distribution must be specified if the epidemic is initialized with fewer subjects than the true population.")
                    } else {
                              if(!(all(names(compartment_dist) %in% epimodel$states) & all(epimodel$states %in% names(compartment_dist)))) {
                                        stop("The names in the compartment_dist vector must match the model compartment names.")
                              }
                    }
          }
          
          if(configs_to_redraw > epimodel$popsize) config_replacement <- TRUE

          # initialize simulation settings list
          epimodel$sim_settings <- list(niter = niter,
                               save_params_every = save_params_every,
                               save_configs_every = save_configs_every,
                               kernel = kernel,
                               cov_mtx = cov_mtx,
                               configs_to_redraw = configs_to_redraw,
                               config_replacement = config_replacement,
                               preferential_sampling = preferential_sampling,
                               init_popsize = init_popsize,
                               compartment_dist = compartment_dist,
                               post_init_params = post_init_params,
                               to_estimation_scale = to_estimation_scale,
                               from_estimation_scale = from_estimation_scale,
                               analytic_eigen = analytic_eigen,
                               ecctmc_method = ecctmc_method,
                               seed = seed,
                               time_limit = time_limit)

          return(epimodel)
}