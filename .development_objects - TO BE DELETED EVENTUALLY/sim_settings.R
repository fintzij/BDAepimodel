# this file contains pieces of code that are used in developing the package but
# are not necessary for the package to work.
init_state = c(S = 45, I = 5, R = 0)
initialization_fcn <- function(){
          init_state <- rmultinom(1, 100, c(0.95, 0.05, 0))
          names(init_state ) <- c("S", "I", "R")
          return(init_state)
}

obstimes = seq(0, 10, by = 1)
meas_vars = "I"

r_meas_process <- function(proc_state, meas_vars, meas_params){
          rbinom(n = 1, size = proc_state$meas_vars, prob = meas_params)
}

epimodel <- init_epimodel(obstimes = 0:10,
                          states = c("S", "I", "R"), 
                          params = c(beta = 0.5, mu = 1, rho = 0.5, p0 = 0.05), 
                          rates = c("beta * I", "mu"), 
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 2), 
                          initialization_fcn = initialization_fcn, 
                          meas_vars = "I",
                          r_meas_process = r_meas_process) 
