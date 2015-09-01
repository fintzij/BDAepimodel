#' Initialize epimodel bookkeeping object.
#' 
#' Initializes a bookkeeping list containing the population level and subject 
#' level bookkeeping objects, measurement, and epidemic process settings. At a 
#' minimum, the epimodel object must be initialized with the states, params, 
#' rates, flow, and either the vector of observation times or a matrix with
#' data.
#' 
#' @param states character vector with names of compartments corresponding 
#'   exactly to the names of compartments used in the \code{*_rates} arguments.
#' @param params character vector of process parameters with element names 
#'   corresponding exactly to names of parameters used in the \code{rates} 
#'   arguments.
#' @param rates list of functions or a character vector with strings specifying 
#'   the flow rates between compartments with compartments corresponding to 
#'   names given by \code{states}.
#' @param flow numeric matrix with reactions as columns and where each row in a 
#'   column indicates the change in the size of the compartment on the subject 
#'   level.
#' @param dat matrix of dimension \code{number of observation times \emph{x} 
#'   number of measured compartments}. \code{dat} must have one column with 
#'   observation times (numeric and strictly increasing), whose the name of 
#'   which is given by the \code{times} variable. The matrix must be sorted in 
#'   ascending time order.
#' @param time_var string indicating the name of the variable coding observation
#'   times in \code{dat}.
#' @param popsize size of the population.
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
#' @param covar optional numeric matrix of time varying covariates (currently 
#'   not working).
#' @param tcovar numeric vector of times for time-varying covariates (currently 
#'   not working).
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
init_epimodel <- function(states, params, rates, flow, dat = NULL, time_var = NULL, obstimes = NULL, popsize = NULL, pop_mat = NULL, subj_mat = NULL, obs_mat = NULL, meas_vars = NULL, r_meas_process = NULL, d_meas_process = NULL, covar = NULL, tcovar = NULL, rprior = NULL, dprior = NULL, to_estimation_scale = NULL, from_estimation_scale = NULL, init_state = NULL) {
          
          
          # user must specify states, parameters, flow, and rates at a minimum. 
          if(missing(states)){
                    stop(sQuote("states"), "must be specified")
          }
          
          if(missing(params)){
                    stop(sQuote("params"), "must be specified")
          }
          
          if(missing(rates)){
                    stop(sQuote("rates"), "must be specified")
          }
          
          if(missing(flow)){
                    stop(sQuote("flow"), "must be specified")
          }
          
          if(missing(dat) & missing(obstimes)){
                    stop("Either the observation times or a dataset must be specified")
          }
          
          # check that the user has indicated which is the times variable if a 
          # dataset has been specified. 
          if(!is.null(dat) & is.null(time_var)){
                    stop(sQuote("time_var"), "must be specified for the dataset")
          }
          
          # check that the user has specified the times along with the covariates
          if(!is.null(covar) & is.null(tcovar)){
                    stop(sQuote("tcovar"), "must be specified along with covariates")
          }
          
          # check that if the one of the transformation arguments was provided,
          # then its inverse should be as well.
          if(is.null(to_estimation_scale) & !is.null(from_estimation_scale)){
                stop(sQuote("to_estimation_scale"),"must be provided if from_estimation_scale is specified")    
          }
          
          if(is.null(from_estimation_scale) & !is.null(to_estimation_scale)){
                    stop(sQuote("from_estimation_scale"),"must be provided if to_estimation_scale is specified")
          }
          
          # if flow rates were provided as a character vector, extract the rates
          # to instatiate a function list
          if(is.character(rates)) {
                    rates <- extract_rate_fcns(rates = rates, states = states, param_names = names(params))
          }
          
          # if a dataset was provided instead of obstimes, instatiate obstimes. 
          # Note that obstimes and dat cannot both be NULL and that time_var is
          # provided if dat is provided.
          if(is.null(obstimes)) {
                    obstimes <- dat[,time_var]
          }
          
          # if a vector of initial compartment counts is provided, deduce the population size
          if(!is.null(init_state)) {
                    popsize = sum(init_state)
          }
          
          #initialize list object
          epimodel <- list(dat = dat,
                           time_var = time_var,
                           obstimes = obstimes,
                           popsize = popsize,
                           states = states,
                           params = params,
                           rates = rates,
                           flow = flow, 
                           pop_mat = pop_mat,
                           subj_mat = subj_mat, 
                           obs_mat = obs_mat,
                           meas_vars = meas_vars,
                           r_meas_process = r_meas_process,
                           d_meas_process = d_meas_process,
                           covar = covar,
                           tcovar = tcovar,
                           rprior = rprior,
                           dprior = dprior,
                           to_estimation_scale = to_estimation_scale,
                           from_estimation_scale = from_estimation_scale)
          
          # if the meas_vars was specified, initialize obs_mat
          if(!is.null(epimodel$meas_vars)){
                    epimodel$obs_mat <- init_obs_mat(epimodel = epimodel)
          }
          
          
          # if the initial vector of compartment counts was provided, instatiate
          # the population and subject level matrices
          if(!is.null(init_state)){
                    epimodel$subj_mat <- init_subj_mat(epimodel = epimodel, init_state = init_state)
                    epimodel$pop_mat <- init_pop_mat(init_state = init_state, tmax = max(obstimes)) 
          }
          
          # return bookkeeping list
          return(epimodel)
}