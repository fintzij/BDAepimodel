library(BDAepimodel)

context("Checking and instatiating a new rate matrix and eigen decomposition")

test_that("It is detected whether a missing rate matrix is needed", {

          set.seed(52787)
          require(BDAepimodel)
          
          popsize <- 10
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
          }
          
          
          epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.5),
                                    popsize = popsize,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = rnorm(1, 0.5, 1e-6), mu = rnorm(1, 1, 1e-6), rho = 0.5,  S0 = 0.5, I0 = 0.5, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          epimodel <- prepare_epimodel(epimodel)
          
          epimodel$pop_mat[6, 4:6] <- c(0, 9, 1)
          
          # remove the rate matrix and decompositions for keys 0 and 3
          epimodel$keys <- retrieveKeys(1:epimodel$ind_final_config, epimodel$irm_key_lookup, epimodel$pop_mat, epimodel$index_state_num)
          
          expect_true(any(epimodel$keys == 0))
          # indices of .keys vector for which keys are missing
          missing_keys <- unique(epimodel$pop_mat[which(epimodel$keys == 0),"I"])
          
          expect_equal(missing_keys, 9)
          
})

test_that("Missing rate matrices and eigen decompositions are instatiated", {
          
          
          set.seed(52787)
          require(BDAepimodel)
          
          popsize <- 10
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
          }
          
          
          epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.5),
                                    popsize = popsize,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = rnorm(1, 0.5, 1e-6), mu = rnorm(1, 1, 1e-6), rho = 0.5,  S0 = 0.5, I0 = 0.5, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          epimodel <- prepare_epimodel(epimodel)
          
          epimodel$pop_mat[6, 4:6] <- c(0, 9, 1)
          
          # remove the rate matrix and decompositions for keys 0 and 3
          epimodel$keys <- retrieveKeys(1:epimodel$ind_final_config, epimodel$irm_key_lookup, epimodel$pop_mat, epimodel$index_state_num)
          
          epimodel <- append_missing(epimodel)
          
          expect_equal(epimodel$irm[,,10], matrix(c(-9 * epimodel$params["beta"], 9 * epimodel$params["beta"], 0, 
                                                   0, -epimodel$params["mu"], epimodel$params["mu"],
                                                   0,0,0), byrow = T, nrow = 3))
          
          expect_equal(epimodel$eigen_vectors[,,10], eigen(epimodel$irm[,,10])$vectors)
})