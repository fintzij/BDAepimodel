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
#'   arguments. The params vector must contain a number of parameters equal to 
#'   the number of states for the categorical distribution for the state of a 
#'   subject at t0. The names of these parameters must end in 0 - e.g. for the 
#'   SIR model: S0, I0, R0. A value of zero for one of the initial distribution 
#'   parameters will be interpretted as a model assumption that subjects cannot 
#'   occupy that state at time t0.
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
#' @param initdist_prior vector of parameters for dirichlet prior for default 
#'   categorical distribution over initial distribution at t0. Defaults to flat 
#'   prior unless otherwise specified.
#' @param meas_vars character vector specifying which compartments are measured.
#' @param incidence_vars character vector specifying which of the measured
#'   compartments contain incidence data rather than prevalence data.
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
#' @param sim_settings bookkeeping list for simulation settings.
#'   
#' @return list containing bookkeeping objects and model configuration objects.
#' @export
#' 
#' @examples epimodel <- init_epimodel(states = c("S", "I", "R"),
#' params = c(beta = 0.02, mu = 1, gamma = 0.5, rho = 0.5, p0 = 0.05),
#' rates = c("beta * I", "mu", "gamma"),
#' flow = matrix(c(-1, 1, 0, 0, -1, 1, 1, 0, -1), ncol = 3, byrow = T))
#' 
init_epimodel <- function(states, params, rates, flow, dat = NULL, time_var = NULL, obstimes = NULL, popsize = NULL, config_mat = NULL, obs_mat = NULL, initdist_prior = NULL, meas_vars = NULL, incidence_vars = NULL, r_meas_process = NULL, d_meas_process = NULL, covar = NULL, tcovar = NULL, sim_settings = NULL) {


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
          
          if(!is.null(incidence_vars) & !all(incidence_vars %in% meas_vars)) {
                    stop(sQuote("The incidence variables must be a subset of the measured variables"))
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

          # check that all of the parameters for the initial distribution have
          # been specified and sum to 1
          if(!all(paste0(states,0) %in% names(params))) {
                    stop("the probability that an individual is in a state at time t0 must be specified for all states.")
          }

          if(sum(params[paste0(states,0)]) != 1) {
                    stop("The vector of probabilities for the state at time t0 must sum to 1.")

          }

          # extract the rate_map
          rate_map <- init_rate_map(rates = rates, states = states)

          # extract the rates to instatiate a function list
          rates <- extract_rate_fcns(rates = rates, states = states, param_names = names(params))

          # if a dataset was provided instead of obstimes, instatiate obstimes.
          # Note that obstimes and dat cannot both be NULL and that time_var is
          # provided if dat is provided.
          if(is.null(obstimes) & !is.null(dat)) {
                    obstimes <- dat[,time_var]
          }

          # identify the parameters for the initial distribution
          initdist_params     <- names(params)[names(params) %in% paste0(states, 0)]

          # get indices of parameters to increment
          initdist_param_inds <- which(params[initdist_params] != 0)


          # set the parameters for the default flat prior for the initial
          # distribution if none provided
          if(is.null(initdist_prior)) {
                    initdist_prior <- rep(1, sum(params[initdist_params] != 0))
          } else {
                    if(length(initdist_prior) != length(initdist_params)) {
                              stop("Dirichlet prior parameters should have length equal to the number of non-zero initial probabilities")
                    }
          }

          # set up the function for sampling from the initial distribution
          # subject level simulation of initial state at time t0
          r_initdist <- function(epimodel) {
                    sample.int(epimodel$num_states, 1, prob = epimodel$params[epimodel$initdist_params])
          }

          # set the function to evaluate the initial distribution for a single subject
          d_initdist <- function(state, params, log = TRUE) {

                    if(log == TRUE) {
                              sum(log(params[state == 1:length(params)]))
                    } else {
                              prod(params^(state == 1:length(params)))
                    }

          }

          # set up the transition kernel for parameters of the initial distribution
          initdist_kernel <- function(epimodel) {

                    params <- epimodel$params
                    
                    # calculate parameters of posterior
                    posterior_params <- epimodel$initdist_prior + epimodel$pop_mat[1, epimodel$states[epimodel$initdist_param_inds]]

                    # sample new parameter values
                    params[epimodel$initdist_params[epimodel$initdist_param_inds]] <- MCMCpack::rdirichlet(1, alpha = posterior_params)
                    
                    return(params)

          }

          # ensure that the flow matrix is an integer matrix
          flow <- matrix(as.integer(flow), nrow = nrow(flow), ncol = ncol(flow))
          colnames(flow) <- states


          #initialize list object
          epimodel <- list(dat = dat,
                           time_var = time_var,
                           obstimes = obstimes,
                           popsize = popsize,
                           states = states,
                           params = params,
                           rates = rates,
                           rate_map = rate_map,
                           flow = flow,
                           config_mat = config_mat,
                           obs_mat = obs_mat,
                           d_initdist = d_initdist,
                           r_initdist = r_initdist,
                           initdist_params = initdist_params,
                           initdist_param_inds = initdist_param_inds,
                           initdist_prior = initdist_prior,
                           initdist_kernel = initdist_kernel,
                           meas_vars = meas_vars,
                           incidence_vars = incidence_vars,
                           r_meas_process = r_meas_process,
                           d_meas_process = d_meas_process,
                           covar = covar,
                           tcovar = tcovar,
                           sim_settings = sim_settings)

          # lookup table for mapping states codes to rates, and vector of states
          # that index the rates
          epimodel$state_lookup <- init_state_lookup(epimodel)
          epimodel$index_states <- colnames(epimodel$rate_map)[apply(epimodel$rate_map, 2, function(x) {1 %in% x})]
          epimodel$index_state_num <- which(epimodel$states %in% epimodel$index_states)

          # if the time_var argument was not supplied, default to "time"
          if(is.null(epimodel$time_var)){
                    epimodel$time_var   <- "time"
          }

          # set the indexing for the observed and measured variables in the
          # observation matrix
          epimodel$meas_vars_aug       <- paste0(epimodel$meas_vars, "_augmented")
          epimodel$meas_vars_obs       <- paste0(epimodel$meas_vars, "_observed")

          # identify the unmeasured compartments, constants for number of
          # measured and unmeasured states
          epimodel$unmeasured_vars     <- setdiff(epimodel$states, epimodel$meas_vars)
          epimodel$num_unmeasured      <- length(epimodel$unmeasured_vars)
          epimodel$num_measured        <- length(epimodel$meas_vars)
          epimodel$num_states          <- length(epimodel$states)

          # set the index matrix for where to insert rates within an irm from the flow matrix
          epimodel$flow_inds <- matrix(nrow = nrow(epimodel$flow), ncol = 2)

          for(j in 1:nrow(epimodel$flow_inds)) {
                  epimodel$flow_inds[j, ] <- c(which(epimodel$flow[j, ] == -1), which(epimodel$flow[j, ] == 1))
          }

          # return bookkeeping list
          return(epimodel)
}