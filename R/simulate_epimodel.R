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
#' @param epimodel list of bookkeeping and model objects.
#' @inheritParams init_epimodel
#' @param return_config option specifying whether to return a matrix of
#'   subject-level trajectories. Defaults to \code{TRUE}.
#' @param trim option specifying whether to trim the simulated data if the
#'   epidemic ends before the final observation time. Defaults to \code{TRUE}.
#' @param lump option specifying whether to simulate the next event directly on
#'   the extended state space of individuals, or to simulate on the lumped state
#'   space of compartment counts, and then choose a subject uniformly at random
#'   for each event from the set of subjects in the risk pool for that event. A
#'   value of TRUE, corresponding to the lumped state space, is appropriate when
#'   there is no subject level heterogeneity and is significantly faster in even
#'   moderate sized populations.
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
#' @examples # Set up initialization funtion and measurement process function
#'
#' r_initdist <- function(params) {
#'                            sample.int(3, 1, prob = c(1-params["rho"], params["rho"], 0))
#'                            }
#'
#' r_meas_process <- function(state, meas_vars, params){
#'                            obs_vec <- rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
#'                            names(obs_vec) <- meas_vars
#'                            return(obs_vec)
#'                            }
#'
#' # initialize epimodel
#' epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.05),
#'                           states = c("S", "I", "R"),
#'                           params = c(beta = 0.02, mu = 1, gamma = 0.5, rho = 0.5, p0 = 0.05),
#'                           rates = c("beta * I", "mu", "gamma"),
#'                           flow = matrix(c(-1, 1, 0, 0, -1, 1, 1, 0, -1), ncol = 3, byrow = T),
#'                           meas_vars = "I",
#'                           r_meas_process = r_meas_process,
#'                           popsize = 200,
#'                           r_initdist = r_initdist)
#'
#' # simulate from model
#' epimodel <- simulate_epimodel(epimodel)
#'
simulate_epimodel <- function(epimodel, obstimes = NULL, init_state = NULL, meas_vars = NULL, r_meas_process = NULL, trim = TRUE, lump = TRUE, ensure_pos_rates = FALSE){

          # check function arguments and issue warnings/errors

          if(is.null(init_state) & (is.null(epimodel$r_initdist) & is.null(epimodel$popsize))) {
                    stop("Either the initial state vector or a function to simulate it must be supplied within the epimodel object along with the population size.")
          }

          if(is.null(epimodel$meas_vars) && is.null(meas_vars)){
                    stop(sQuote("meas_vars"), "must be specified either within the epimodel list or as an argument to the simulation function.")
          }

          if(is.null(init_state) & is.null(epimodel$popsize)) {
                    stop("If an initial state vector is not provided, the population size must be specified.")
          }

          if(is.null(epimodel$r_meas_process) && is.null(r_meas_process)){
                    stop(sQuote("r_meas_process"), "must be specified either within the epimodel list or as an argument to the simulation function.")
          }

          if(is.null(epimodel$obstimes) && is.null(obstimes)) {
                    stop("The observation times must be supplied, either within the epimodel list or as an argument to the simulation function.")
          }
          
          # generate initial state vector if initial distribution is not supplied
          if(is.null(init_state)) {
                    init_config <- replicate(n = epimodel$popsize, expr = epimodel$r_initdist(epimodel))
                    init_state <- numeric(length = length(epimodel$states)); names(init_state) <- epimodel$states
                    for(k in 1:length(init_state)) {
                              init_state[k] <- sum(init_config == k)
                    }
          } else {
                  init_config <- rep.int(1:length(epimodel$states), times = init_state)
          }

          # get population size
          epimodel$popsize <- sum(init_state)

          # initialize bookkeeping matrix
          if(is.null(epimodel$obstimes) && !is.null(obstimes)){
                    epimodel$obstimes <- obstimes
          }

          if(is.null(epimodel$meas_vars) & !is.null(meas_vars)) {
                  epimodel$meas_vars <- meas_vars
          }

          # initialize the configuration and population matrices
          epimodel <- init_config_mat(epimodel, init_state = init_state, t0 = min(epimodel$obstimes), tmax = max(epimodel$obstimes))

          # initialize simulation
          subj_config <- epimodel$config_mat[1, , drop = FALSE]
          pop_config  <- epimodel$pop_mat[1, , drop = FALSE]

          if(lump == TRUE) {
                    rate_mat <- matrix(1, nrow = 1, ncol = nrow(epimodel$flow)) # initialize rate matrix
                    dt_mat <- matrix(Inf, nrow = 1, ncol = nrow(epimodel$flow)) # initialize dt matrix

          } else if(lump == FALSE) {
                    rate_mat <- matrix(1, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize rate matrix
                    dt_mat <- matrix(Inf, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize dt matrix
          }

          config_list <- list(pop_config = pop_config, subj_config = subj_config, rate_mat = rate_mat, dt_mat = dt_mat, keep_going = TRUE)

          # control objects for the simulation
          tmax <- max(epimodel$obstimes)

          # expand the pop_mat and initialize row index
          final_row <- epimodel$pop_mat[2, ]
          epimodel$pop_mat <- rbind(epimodel$pop_mat[1,],
                                    matrix(0.0,
                                           nrow = epimodel$popsize * nrow(epimodel$flow),
                                           ncol = ncol(epimodel$pop_mat)
                                    ))
          row_ind <- 1

          # simulate until tmax or no more events
          while((config_list$pop_config[,1,drop=F] <= tmax) && config_list$keep_going) {

                    config_list <- sim_one_event(config_list, epimodel, lump)
                    row_ind <- row_ind + 1

                    if(row_ind > nrow(epimodel$pop_mat)) {
                              epimodel$pop_mat <- rbind(epimodel$pop_mat,
                                                        matrix(0.0,
                                                               nrow = epimodel$popsize*nrow(epimodel$flow),
                                                               ncol = ncol(epimodel$pop_mat)))
                    }

                    if(config_list$keep_going){# insert the new row

                              epimodel$pop_mat[row_ind,] <- config_list$pop_config

                    } else { # replace the final row
                            epimodel$pop_mat[row_ind,] <- config_list$pop_config
                            # epimodel$config_mat[nrow(epimodel$config_mat),] <- config_list$subj_config
                    }
          }

          epimodel$pop_mat <- epimodel$pop_mat[1:row_ind,]

          # initialize obs_mat
          epimodel$obs_mat <- init_obs_mat(epimodel)

          # sample from the measurement process
          epimodel$obs_mat[, paste(epimodel$meas_vars, "_observed", sep = "")] <-
                    epimodel$r_meas_process(epimodel$obs_mat,
                                            paste(epimodel$meas_vars, "_augmented", sep = ""),
                                            epimodel$params)

          # instatiate data matrix
          epimodel$dat <- epimodel$obs_mat[, c("time", paste(epimodel$meas_vars, "_observed", sep = ""))]
          colnames(epimodel$dat) <- c(epimodel$time_var, epimodel$meas_vars)

          # if trim is TRUE, then trim off the excess 0 measurements in the case
          # when the epidemic died off before the final observation time.
          if((trim == TRUE) && all(config_list$rate_mat == 0)) {
                    # find the index of the first observation time where the rates are 0
                    end_ind <- which(epimodel$obstimes > epimodel$pop_mat[nrow(epimodel$pop_mat) - 1, "time"])[1]

                    # trim the observation times, observation matrix, and data matrix
                    epimodel$obstimes <- epimodel$obstimes[1:end_ind]
                    epimodel$obs_mat <- epimodel$obs_mat[1:end_ind, ]
                    epimodel$dat <- epimodel$dat[1:end_ind, ]
          }

          epimodel$init_config <- epimodel$config_mat[1,]
          epimodel$config_mat <- NULL

          # if(return_config == FALSE) {
          #           epimodel$pop_mat <- NULL
          # }

          return(epimodel)

}