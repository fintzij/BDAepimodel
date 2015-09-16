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
#' @param params numeric vector of process parameters with element names 
#'   corresponding exactly to names of parameters used in the \code{rates} 
#'   arguments.
#' @param rates character vector with strings specifying the flow rates between 
#'   compartments. Each variable in a rate must correspond to the name of one of
#'   the \code{states}, \code{params}, or \code{covars} (covariates not yet 
#'   implemented).
#' @param flow numeric matrix of dimension \code{number of transitions \emph{x}
#'   number of compartments}. Each row corresponds to a possible transition, and
#'   each column in a row has element 1 to indicate an entry to that
#'   compartment, -1 to indicate an exit, and 0 for no change in the size of the
#'   compartment on the subject level. 
#' @param dat matrix of dimension \code{number of observation times \emph{x} 
#'   number of measured compartments}. \code{dat} must have one column with 
#'   observation times (numeric and strictly increasing), whose the name of 
#'   which is given by the \code{times} variable. The matrix must be sorted in 
#'   ascending time order.
#' @param time_var string indicating the name of the variable coding observation
#'   times in \code{dat}.
#' @param popsize size of the population.
#' @param config_mat population-level bookkeeping matrix for compartment counts 
#'   and the configuration of individuals at times of state transition.
#' @param obs_mat population-level bookkeeping matrix for compartment counts at 
#'   observation times.
#' @param r_initdist function to simulate the subject-level initial distribution
#'   at time t0.
#' @param d_initdist function to evaluate the subject-level distribution at time
#'   t0.
#' @param meas_vars character vector specifying which compartments are measured.
#' @param r_meas_process function to simulate from the measurement process, with
#'   named arguments \code{state}, \code{meas_vars}, and \code{params}, which 
#'   are the current state of the process (given as a named character vector 
#'   with element names corresponding exactly to the names of compartments), a 
#'   vector with the names of the states to be measured, and a named vector of 
#'   process parameters. The function should return noisy measurements of the 
#'   process in a vector with named elements corresponding to each of the 
#'   measured compartments specified in \code{meas_vars}. The function must 
#'   index into the state vector and parameter vector using element names.
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
#' @examples epimodel <- init_epimodel(states = c("S", "I", "R"), 
#' params = c(beta = 0.02, mu = 1, gamma = 0.5, rho = 0.5, p0 = 0.05),
#' rates = c("beta * I", "mu", "gamma"),
#' flow = matrix(c(-1, 1, 0, 0, -1, 1, 1, 0, -1), ncol = 3, byrow = T))
#' 
init_epimodel <- function(states, params, rates, flow, dat = NULL, time_var = NULL, obstimes = NULL, popsize = NULL, config_mat = NULL, obs_mat = NULL, r_initdist = NULL, d_initdist = NULL, meas_vars = NULL, r_meas_process = NULL, d_meas_process = NULL, covar = NULL, tcovar = NULL, rprior = NULL, dprior = NULL, to_estimation_scale = NULL, from_estimation_scale = NULL) {
          
          
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
          
          if(ncol(flow) != length(states)){
                    stop(sQuote("flow"), "must have number of columns equal to the number of states")
          } else{
                    colnames(flow) <- states
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
          
          # extract the rates to instatiate a function list
          rates <- extract_rate_fcns(rates = rates, states = states, param_names = names(params))
          
          # if a dataset was provided instead of obstimes, instatiate obstimes. 
          # Note that obstimes and dat cannot both be NULL and that time_var is
          # provided if dat is provided.
          if(is.null(obstimes) & !is.null(dat)) {
                    obstimes <- dat[,time_var]
          }
          
                    
          #initialize list object
          epimodel <- structure(list(dat = dat,
                           time_var = time_var,
                           obstimes = obstimes,
                           popsize = popsize,
                           states = states,
                           state_lookup = init_state_lookup(states),
                           params = params,
                           rates = rates,
                           flow = flow, 
                           config_mat = config_mat,
                           obs_mat = obs_mat,
                           r_initdist = r_initdist,
                           d_initdist = d_initdist,
                           meas_vars = meas_vars,
                           r_meas_process = r_meas_process,
                           d_meas_process = d_meas_process,
                           covar = covar,
                           tcovar = tcovar,
                           rprior = rprior,
                           dprior = dprior,
                           to_estimation_scale = to_estimation_scale,
                           from_estimation_scale = from_estimation_scale), class = "epimodel")
          
          # if the time_var argument was not supplied, default to "time"
          if(is.null(epimodel$time_var)){
                    epimodel$time_var <- "time"
          }
          
          
          # return bookkeeping list
          return(epimodel)
}