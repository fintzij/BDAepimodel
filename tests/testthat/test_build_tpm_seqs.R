library(BDAepimodel)

context("Constructing the tpms and tpm product sequences")

test_that("tpms are correctly computed", {
          
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
                                    params =  c(beta = 0.9, mu = 1, rho = 0.5,  S0 = 0.5, I0 = 0.5, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          epimodel <- prepare_epimodel(epimodel)
          
          tpmSeqs(tpms = epimodel$tpms, pop_mat = epimodel$pop_mat, eigen_vals = epimodel$eigen_values, eigen_vecs = epimodel$eigen_vectors, inverse_vecs = epimodel$inv_eigen_vectors, irm_keys = epimodel$keys)
          
          tpm1 <- eigen(epimodel$irm[,,1])$vectors %*% diag(exp(eigen(epimodel$irm[,,1])$values * (epimodel$pop_mat[2, "time"] - epimodel$pop_mat[1, "time"]))) %*% solve(eigen(epimodel$irm[,,1])$vectors)
          
          tpm6 <- eigen(epimodel$irm[,,3])$vectors %*% diag(exp(eigen(epimodel$irm[,,3])$values * (epimodel$pop_mat[7, "time"] - epimodel$pop_mat[6, "time"]))) %*% solve(eigen(epimodel$irm[,,3])$vectors)
          
          tpm7 <- eigen(epimodel$irm[,,3])$vectors %*% diag(exp(eigen(epimodel$irm[,,3])$values * (epimodel$pop_mat[8, "time"] - epimodel$pop_mat[7, "time"]))) %*% solve(eigen(epimodel$irm[,,3])$vectors)
          
          expect_equal(epimodel$tpms[,,1], tpm1)
          expect_equal(epimodel$tpms[,,6], tpm6)
          expect_equal(epimodel$tpms[,,7], tpm7)
}) 


test_that("tpm product sequences are computed properly", {
          
          
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
                                    params =  c(beta = 0.9, mu = 1, rho = 0.5,  S0 = 0.5, I0 = 0.5, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          epimodel <- prepare_epimodel(epimodel)
          
          tpmSeqs(tpms = epimodel$tpms, pop_mat = epimodel$pop_mat, eigen_vals = epimodel$eigen_values, eigen_vecs = epimodel$eigen_vectors, inverse_vecs = epimodel$inv_eigen_vectors, irm_keys = epimodel$keys)
          
          tpmProdSeqs(tpm_prods = epimodel$tpm_products, tpms = epimodel$tpms, obs_time_inds = epimodel$obs_time_inds)
          
          tpm_prod9 <- epimodel$tpms[,,9] %*% epimodel$tpms[,,10]
          
          expect_equal(epimodel$tpm_products[,,9], tpm_prod9)
          
})
