library(BDAepimodel) 

context("Generating keys to index into the .irm and .eigen environments")

test_that("a subset of keys is generated properly", {
          
          set.seed(52787)
          
          popsize <- 5
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
          }
          
          epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.5),
                                    popsize = popsize,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = 0.9, mu = 1, rho = 0.5, S0 = 0.5, I0 = 0.5, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          epimodel <- prepare_epimodel(epimodel)
          
          expect_equal(as.numeric(retrieveKeys(inds = 5, irm_lookup = epimodel$irm_key_lookup, pop_mat = epimodel$pop_mat, epimodel$index_state_num)), 2)
          expect_equal(as.numeric(retrieveKeys(inds = 5:7, irm_lookup = epimodel$irm_key_lookup, pop_mat = epimodel$pop_mat, epimodel$index_state_num)), c(2, 3, 3))
          expect_equal(as.numeric(retrieveKeys(inds = 1:epimodel$ind_final_config, irm_lookup = epimodel$irm_key_lookup, pop_mat = epimodel$pop_mat, epimodel$index_state_num)), epimodel$irm_key_lookup[match(epimodel$pop_mat[1:epimodel$ind_final_config,"I"], epimodel$irm_key_lookup[,"I"]),"key"])
})
