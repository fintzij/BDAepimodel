# this file contains pieces of code that are used in developing the package but
# are not necessary for the package to work.
popsize <- 200

initialization_fcn <- function(){
          init_state <- as.numeric(rmultinom(1, popsize, c(0.95, 0.05, 0)))
          names(init_state ) <- c("S", "I", "R")
          return(init_state)
}

r_meas_process <- function(state, meas_vars, params){
          rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
}

# R0 = 4, mu = 1, rho = 0.5, p0 = 0.05
epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.5),
                          states = c("S", "I", "R"), 
                          params = c(beta = 0.02, mu = 1, rho = 0.5, p0 = 0.05), 
                          rates = c("beta * I", "mu"), 
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                          meas_vars = "I",
                          r_meas_process = r_meas_process)

init_state <- initialization_fcn()
if(init_state["I"] == 0){
          while(init_state["I"] == 0){
                    init_state <- initialization_fcn()
          }
}

epimodel <- simulate_epimodel(epimodel = epimodel, init_state = init_state, lump = TRUE, trim = FALSE)
