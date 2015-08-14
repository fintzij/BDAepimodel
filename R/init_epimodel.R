#' Initialize epimodel bookkeeping object.
#' 
#' Initializes a bookkeeping list containing the population level and subject 
#' level bookkeeping objects, measurement, and epidemic process settings. At a 
#' minimum, the epimodel object must be initialized with character vectors of 
#' compartment names and parameter names.
#' 
#' @param dat matrix of dimension \code{number of observation times \emph{x} 
#'   number of measured compartments}. \code{dat} must have one column with 
#'   observation times (numeric and strictly increasing), whose the name of 
#'   which is given by the \code{times} variable. The matrix must be sorted in 
#'   ascending time order.
#' @param obstimes string indicating the name of the variable coding observation
#'   times in \code{dat}.
#' @param t_0 time at which the epidemic begins, no later than the first 
#'   observation time.
#' @param popsize size of the population.
#' @param states character vector with names of compartments corresponding 
#'   exactly to the names of compartments used in the \code{*_rates} arguments.
#' @param params named vector of process parameters with element names 
#'   corresponding exactly to names of parameters used in the \code{*_rates} 
#'   arguments.
#' @param pop_mat population-level bookkeeping matrix for compartment counts at 
#'   times of state transition.
#' @param subj_mat subject-level bookkeeping matrix for individual trajectories.
#' @param obs_mat population-level bookkeeping matrix for compartment counts at 
#'   observation times.
#' @param meas_vars character vector specifying which compartments are measured.
#' @param r_meas_process function to simulate from the measurement process, with
#'   named arguments \code{proc_state} and \code{meas_params}, which are the 
#'   current state of the process (given as a named character vector with 
#'   element names corresponding exactly to the names of compartments), and 
#'   \code{meas_params}, a vector of named measurement process parameters. The 
#'   function should return noisy measurements of the process in a vector with 
#'   named elements corresponding to each of the measured compartments specified
#'   in \code{meas_vars}.
#' @param d_meas_process function to evaluate the density of the measurement 
#'   process with arguments and output specified as in \code{r_meas_process}.
#' @param pop_rates list of functions for computing the population level flow 
#'   rates between compartments with compartments corresponding to names given 
#'   by \code{states}.
#' @param pop_flow numeric matrix with reactions as columns and where each row 
#'   in a column indicates the change in the size of the compartment on the 
#'   population level.
#' @param subj_rates list of functions for computing the flow rates between 
#'   compartments with compartments corresponding to names given by 
#'   \code{states}.
#' @param subj_flow numeric matrix with reactions as columns and where each row 
#'   in a column indicates the change in the size of the compartment on the 
#'   subject level.
#' @param covar optional numeric matrix of time varying covariates (currently 
#'   not working).
#' @param tcovar numeric vector of times for time-varying covariates.
#' @param rprior list of functions for simulating from the prior distributions 
#'   of model parameters. List element names should correspond to exactly to the
#'   parameter names given in \code{params}.
#' @param dprior list of functions for evaluating the prior densities of model 
#'   parameters. List element names should correspond to exactly to the 
#'   parameter names given in \code{params}.
#' @param to_estimation_scale list of functions for transforming model 
#'   parameters to the scale on which new parameters should be proposed. List 
#'   element names should correspond to exactly to the parameter names given in 
#'   \code{params}.
#' @param from_estimation_scale list of functions for transforming model 
#'   parameters from the a scale on which new parameters are proposed to the 
#'   scale used to evaluate the process likelihood. List element names should 
#'   correspond to exactly to the parameter names given in \code{params}.
#'   
#' @return list containing bookkeeping objects and model configuration objects.
#' @export
#' 
init_epimodel <- function(dat = NULL, times = NULL, t_0 = 0, popsize = NULL, states, params, pop_mat = NULL, subj_mat = NULL, obs_mat = NULL, meas_vars = NULL, r_meas_process = NULL, d_meas_process = NULL, pop_rates = NULL, pop_flow= NULL, subj_rates = NULL, subj_flow = NULL, covar = NULL, tcovar = NULL, rprior = NULL, dprior = NULL, to_estimation_scale = NULL, from_estimation_scale = NULL) {
          
          # user must specify states and parameters at a minimum. 
          if(missing(states)){
                    stop(sQuote("states"), "must be specified")
          }
          
          if(missing(params)){
                    stop(sQuote("params"), "must be specified")
          }
          
          # check that the user has indicated which is the times variable if a
          # dataset has been specified. Also check this for the time varying covariates
          if(!is.null(dat) & is.null(times)){
                    stop(sQuote("times"), "must be specified")
          }
          
          if(!is.null(covar) & is.null(tcovar)){
                    stop(sQuote("tcovar"), "must be specified")
          }
          
          # check that if the one of the transformation arguments was provided,
          # then its inverse should be as well.
          if(is.null(to_estimation_scale) & !is.null(from_estimation_scale)){
                stop(sQuote("to_estimation_scale"),"must be provided if from_estimation_scale is specified")    
          }
          
          if(is.null(from_estimation_scale) & !is.null(to_estimation_scale)){
                    stop(sQuote("from_estimation_scale"),"must be provided if to_estimation_scale is specified")
          }
          
          #initialize list object
          epimod_mats <- list(dat = dat, 
                              obstimes = obstimes, 
                              t_0 = t_0,
                              popsize = popsize,
                              states = states,
                              params = params,
                              pop_mat = pop_mat, 
                              subj_mat = subj_mat,
                              obs_mat = obs_mat, 
                              meas_vars = meas_vars, 
                              r_meas_process = r_meas_process, 
                              d_meas_process = d_meas_process, 
                              pop_rates = pop_rates,
                              pop_flow = pop_flow, 
                              subj_rates = subj_rates, 
                              subj_flow = subj_flow,
                              covar = covar, 
                              tcovar = tcovar, 
                              rprior = rprior, 
                              dprior = dprior, 
                              to_estimation_scale = to_estimation_scale, 
                              from_estimation_scale = from_estimation_scale)
          
          
          # if a dataset is provided, initialize obs_mat, pop_mat, and subj_mat 
          # n.b. reconsider whether to instatiate subj_mat here or in the
          # simulate_epimodel function
          if(!is.null(dat)){
                    epimodel$obs_mat    <- .init_obs_mat(epimodel)
                    epimodel$pop_mat    <- .init_pop_mat(epimodel)
                    epimodel$subj_mat   <- .init_subj_mat(epimodel) 
          } 
          
          
          # return bookkeeping list
          return(epimod_mats)
}