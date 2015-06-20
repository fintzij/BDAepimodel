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
#' 'GillespieSSA' package so as to facilitate limit confusion for users familiar
#' with that package. 
#' 
#' @param proc_params named vector of process parameters with elements 
#'        corresponding exactly to names of parameters used in the \code{rates} 
#'        vector.
#' @param init_state numeric vector of initial states of compartments the system
#'        with named elements corresponding exactly to the names of compartments
#'        used in the \code{rates} argument.  
#' @param rates character vector of rates of flow rates between compartments
#'        with compartments corresponding to names of elements of 
#'        \code{init_state}.
#' @param flow numeric matrix with reactions as columns and where each row in a
#'        column indicates the change in the size of the compartment. 
#' @param obstimes numeric vector of observation times with first element the 
#'        time at which the process is initialized, and last element the final 
#'        time of observation. The user may specify \code{Inf} as an option, in 
#'        which case the epidemic is completely observed. 
#' @param meas_process function specifying the sampling mechanism of the 
#'        measurement process (e.g. binomial), with named arguments proc_state 
#'        and meas_params, which are the current state of the process (of the 
#'        same form as init_state) and the vector of named measurement process 
#'        parameters. The function should return a noisy measurement.
#' @param meas_params named vector of measurement process parameters with 
#'        elements corresponding exactly to the arguments to \code{meas_process}
#' @param return_trajectories option specifying whether to return a matrix of 
#'        subject-level trajectories. Defaults to \code{FALSE}.
#' @param tmax time until which to simulate the process, possibly infinity.
#' 
#' @return Data frame with columns for the observation times, and draws
#'        from the measurement process at observation times. Optionally, also 
#'        return a data frame with columns for event times and states of the 
#'        epidemic compartments at each event time. 
#'
#' @references {Gillespie DT (1976). “A general method for numerically simulating 
#'        the stochastic time evolution of coupled chemical reactions.” Journal 
#'        of Computational Physics, *22*(4), pp. 403-434.
#' 
#'        Gillespie DT (1977). “Exact stochastic simulation of coupled chemical 
#'        reactions.” The Journal of Physical Chemistry, *81*(25), pp. 2340-
#'        2361.}
#'        
#' 
#' @export
#' 
simulate_epimodel <- function(proc_params, init_state, rates, flow, obstimes, meas_process, meas_params, tmax, return_trajectories = FALSE){
          # extract compartment names
          .compartments <- names(init_state)
          
          # extract parameter names
          .params <- names(proc_params)
          
          # define propensity functions
          .rate_fcns <- extract_rate_fcns(rates = rates, compartments = .compartments, params = .params)
          
          
          ##################################################################################
          ### Completed: function to extract rate functions. 
          ### Next: function to take list of rate functions and current params and state
          ###       in order to compute current rates
          
          # extract process parameters 
          # names
          param_names <- names(proc_params)
          
          # parameters, assigned to names
          for(s in seq_along(proc_params)){
                    assign(param_names[s], proc_params)
          }
          
          # set up data object
          .dat <- init_data(obstimes = obstimes, compartments = compartments)
          
          # set up the population trajectory matrix - i.e. population level count matrix
          pop_mat <- init_pop_traj(compartments = compartments, init_state = init_state)
          
}