library(BDAepimodel)

context("Instatiating the observation matrix")

test_that("The compartment counts at observation times are correct", {
          # initialize epimodel
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          epimodel <- init_epimodel(obstimes = 0:10,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = 0.5, mu = 1, rho = 0.5, S0 = 0.95, I0 = 0.05, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process)
          
          # simulate the epidemic
          set.seed(52787)
          epimodel <- simulate_epimodel(epimodel, init_state = c(S = 7, I = 3, R = 0), lump = TRUE) 
          
          expect_equal(as.numeric(epimodel$obs_mat[,"I_augmented"]), c(3, 2, 5, 3, 1, 1, 0)) 
})
