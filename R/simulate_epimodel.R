#' Simulates a stochastic epidemic model. 
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
#'        time of observation.
#' @param meas_process function specifying the sampling mechanism of the 
#'        measurement process (e.g. binomial), that takes as arguments a named
#'        vector of the same form as \code{init_state}, and the vector of named
#'        measurement process parameters, \code{meas_params}. 
#' @param meas_params named vector of measurement process parameters with 
#'        elements corresponding exactly to the arguments to \code{meas_process}
#'        .