#' Simulates a stochastic epidemic model.
#' 
#' Simulate a realization from a stochastic epidemic model along with a dataset 
#' resulting from a specified measurement process. The (latent) population 
#' trajectory for the epidemic is simulated using Gillespie's first reaction
#' method (Gillespie 1976, 1977), with the trajectory optionally being returned
#' in the output. The dataset that is returned consists of observations of the
#' epidemic that are distributed according to a measurement process, provided as
#' a function by the user. The dynamics of the system are characterized by a
#' named vector of flow rates between compartments and a numeric matrix that
#' specifies the changes to each compartment due to each type of event. The user
#' also specifies a named vector of process parameters, a named vector of
#' initial states, and a set of times at which the epidemic process is sampled.
#' The syntax for calling the function is deliberately similar to that used in
#' the 'GillespieSSA' package so as to facilitate ease of use for users familiar
#' with that package.
#' 
#' @param epimodel list of bookkeeping and model objects. If not supplied, the 
#'   list will be generated.
#' @inheritParams init_epimodel
#' @param initialization_function optional function requiring no inputs that can
#'   be provided in place of init_state to simulate the initial state of the 
#'   system. The function must output a named vector of the same form as 
#'   \code{init_state}.
#' @param return_config option specifying whether to return a matrix of 
#'   subject-level trajectories. Defaults to \code{TRUE}.
#' @param trim option specifying whether to trim the simulated data if the 
#'   epidemic ends before the final observation time. Defaults to \code{TRUE}.
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
simulate_epimodel <- function(epimodel, init_state = NULL, initialization_fcn = NULL, return_config = TRUE, trim = TRUE){
          
          # check function arguments and issue warnings/errors
          
          if(is.null(init_state) & is.null(initialization_fcn) & is.null(epimodel$init_state)) {
                    stop("Either the initial state vector or a function to simulate it must be supplied.")
          }
          
          if(!is.null(init_state) & !is.null(initialization_fcn)) {
                    stop("Only one of the initial state vector and an initialization function may be specified.")
          }
          
          if(is.null(epimodel$meas_vars)){
                    stop(sQuote("meas_vars"), "must be specified within the epimodel list.")
          }
          
          if(is.null(epimodel$r_meas_process)){
                    stop(sQuote("r_meas_process"), "must be specified within the epimodel list.")
          }
          
          # generate initial state vector if initialization function is supplied
          if(!is.null(initialization_fcn)){
                    init_state <- initialization_fcn()
          }
          
          # get population size 
          epimodel$popsize <- sum(init_state)
          
          # initialize bookkeeping matrix
          epimodel$config_mat <- init_config_mat(init_state = init_state, t0 = min(epimodel$obstimes), tmax = max(epimodel$obstimes))

          # initialize simulation
          config <- epimodel$config_mat[epimodel$config_mat[,"time"] == min(epimodel$obstimes), , drop = FALSE]
          rate_mat <- matrix(1, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize rate matrix
          dt_mat <- matrix(Inf, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize dt matrix 
          
          config_list <- list(config = config, rate_mat = rate_mat, dt_mat = dt_mat, keep_going = TRUE)

          # simulate until tmax or no more events
          while(config_list$keep_going) {
                    if(config_list$config[,"I"] == 0) break
                    config_list <- sim_one_event(config_list, epimodel)
                    
                    if(config_list$keep_going){# insert the new row

                              epimodel$config_mat <-insert_row(mat = epimodel$config_mat, 
                                                               row_num = nrow(epimodel$config_mat),
                                                               vec = config_list$config)
                    } else { # replace the final row
                              epimodel$config_mat[nrow(epimodel$config_mat),] <- config_list$config
                    }
          }
          
          # initialize obs_mat
          epimodel$obs_mat <- init_obs_mat(epimodel)
          
          # sample from the measurement process
          epimodel$obs_mat[, paste(epimodel$meas_vars, "_observed", sep = "")] <- 
                    r_meas_process(epimodel$obs_mat, paste(epimodel$meas_vars, "_truth", sep = ""), epimodel$params)
          
          # instatiate data matrix
          epimodel$dat <- epimodel$obs_mat[, c("time", paste(epimodel$meas_vars, "_observed", sep = ""))]
          colnames(epimodel$dat) <- c(epimodel$time_var, epimodel$meas_vars)
          
          # if trim is TRUE, then trim off the excess 0 measurements in the case
          # when the epidemic died off before the final observation time.
          if((trim == TRUE) && all(config_list$rate_mat == 0)) {
                    # find the index of the first observation time where the rates are 0
                    end_ind <- which(epimodel$obstimes > epimodel$config_mat[nrow(epimodel$config_mat) - 1, "time"])[1]
                    
                    # trim the observation times, observation matrix, and data matrix
                    epimodel$obstimes <- epimodel$obstimes[1:end_ind]
                    epimodel$obs_mat <- epimodel$obs_mat[1:end_ind, ]
                    epimodel$dat <- epimodel$dat[1:end_ind, ]
                    
                    # reset the final observation time in the configuration matrix
                    epimodel$config_mat[nrow(epimodel$config_mat),"time"] <- max(epimodel$obstimes)
          }
          
          if(return_config == FALSE) {
                    epimodel$config_mat <- NULL
          }
          
          return(epimodel)

}