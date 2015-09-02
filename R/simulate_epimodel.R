#' Simulates a stochastic epidemic model.
#' 
#' Simulate a realization from a stochastic epidemic model along with a dataset 
#' resulting from a specified measurement process. The (latent) population 
#' trajectory for the epidemic is simulated using Gillespie's direct method 
#' (Gillespie 1976, 1977), with the trajectory optionally being returned in the 
#' output. The dataset that is returned consists of observations of the epidemic
#' that are distributed according to a measurement process, provided as a 
#' function by the user. The dynamics of the system are characterized by a named
#' vector of flow rates between compartments and a numeric matrix that specifies
#' the changes to each compartment due to each type of event. The user also 
#' specifies a named vector of process parameters, a named vector of initial 
#' states, and a set of times at which the epidemic process is sampled. The 
#' syntax for calling the function is deliberately similar to that used in the 
#' 'GillespieSSA' package so as to facilitate ease of use for users familiar 
#' with that package.
#' 
#' @param epimodel list of bookkeeping and model objects. If not supplied, the 
#'   list will be generated.
#' @inheritParams init_epimodel
#' @param return_trajecs option specifying whether to return a matrix of 
#'   subject-level trajectories. Defaults to \code{TRUE}.
#' @param trim option specifying whether to trim the simulated data if the epidemic ends before the final observation time. Defaults to \code{TRUE}.
#'   
#' @return Data frame with columns for the observation times, and draws from the
#'   measurement process at observation times. Optionally, also return a data 
#'   frame with columns for event times and states of the epidemic compartments 
#'   at each event time.
#'   
#' @references {Gillespie DT (1976). “A general method for numerically 
#'   simulating the stochastic time evolution of coupled chemical reactions.” 
#'   Journal of Computational Physics, *22*(4), pp. 403-434.
#'   
#'   Gillespie DT (1977). “Exact stochastic simulation of coupled chemical 
#'   reactions.” The Journal of Physical Chemistry, *81*(25), pp. 2340- 2361.}
#'   
#' @export
#' 
simulate_epimodel <- function(epimodel, init_state = NULL, initialization_fcn = NULL, return_trajecs = TRUE, trim = TRUE){
          
          # check function arguments and issue warnings/errors
          
          if(is.null(init_state) & is.null(initialization_fcn) & is.null(epimodel)) {
                    stop("Either the initial state vector or a function to simulate it must be supplied.")
          }
          
          if(!is.null(init_state) & !is.null(initialization_fcn)) {
                    stop("Only one of the initial state vector and an initialization function may be specified.")
          }
          
          if(!is.null(epimodel) & is.null(epimodel$meas_vars)){
                    stop(sQuote("meas_vars"), "must be specified within the epimodel list.")
          }
          
          if(!is.null(epimodel) & is.null(epimodel$r_meas_process)){
                    stop(sQuote("r_meas_process"), "must be specified within the epimodel list.")
          }
          
          # generate initial state vector if initialization function is supplied
          if(!is.null(initialization_fcn)){
                    init_state <- initialization_fcn()
          }
          
          # extract compartment and process parameter names
          param_names         <- names(params)
          states              <- names(init_state)

          
          # compute auxilliary quantities
          popsize             <- sum(init_state)
          
          
          # initialize bookkeeping object if one is not supplied. If an epimodel list is supplied, it must contain at a minimum the initialized subject, population, and observation matrices, along with states, parameters, flow, 
          if(missing(epimodel)){
                    epimodel <- init_epimodel(obstimes = obstimes, states = states, params = params, rates = rates, flow = flow, init_state = init_state, r_meas_process = r_meas_process)
          }
          
          
          # initialize simulation
          t <- min(epimodel$obstimes)
          status_vec <- rep(states, init_state)
          subj_hazards <- 
          
          return(epimodel)

}