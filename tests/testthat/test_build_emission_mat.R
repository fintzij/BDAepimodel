library(BDAepimodel)

context("Constructing the matrix of emission probabilities")

test_that("Emission probabilities are computed correctly", {
          set.seed(52787)
          
          popsize <- 5
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
          }
          
         
          # R0 = 4, mu = 1, rho = 0.5, p0 = 0.05
          epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.5),
                                    popsize = popsize,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = 1, mu = 1, rho = 0.5, S0 = 0.5, I0 = 0.5, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          .epimodel <- list2env(epimodel, parent = emptyenv(), hash = TRUE)
          
          .epimodel$unmeasured_vars     <- setdiff(.epimodel$states, .epimodel$meas_vars)
          .epimodel$num_unmeasured      <- length(.epimodel$unmeasured_vars)
          .epimodel$num_measured        <- length(.epimodel$meas_vars)
          
          .epimodel$nobs                <- length(.epimodel$obstimes)
          .epimodel$.emission_mat       <- matrix(0, nrow = length(.epimodel$states), ncol = length(.epimodel$obstimes), dimnames = list(.epimodel$states, .epimodel$obstimes))
          
          build_emission_mat(.epimodel)
          
          expect_equal(unname(.epimodel$.emission_mat["S",]), dbinom(.epimodel$obs_mat[,"I_observed"], .epimodel$obs_mat[,"I_augmented"], 0.5))
          expect_equal(unname(.epimodel$.emission_mat["R",]), dbinom(.epimodel$obs_mat[,"I_observed"], .epimodel$obs_mat[,"I_augmented"], 0.5))
          expect_equal(unname(.epimodel$.emission_mat["I",]), dbinom(.epimodel$obs_mat[,"I_observed"], .epimodel$obs_mat[,"I_augmented"] + 1, 0.5))
          
})