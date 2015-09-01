# this file contains pieces of code that are used in developing the package but
# are not necessary for the package to work.
init_state = c(S = 45, I = 5, R = 0)

obstimes = seq(0, 10, by = 1)
meas_vars = "I"

meas_process <- function(proc_state, meas_params){
          rbinom(n = 1, size = proc_state$I, prob = meas_params)
}

init_epimodel(obstimes = 0:10,states = c("S", "I", "R"), params = c(beta = 0.5, mu = 1, rho = 0.5, p0 = 0.05), rates = c("beta * S * I", "mu * I"), flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 2), init_state = c(S = 40, I = 5, R = 5)) 
