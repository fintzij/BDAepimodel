library(BDAepimodel)

context("Test that re-inserting a path reproduces the original configuration matrix")

test_that("The original configuration matrix is restored when a path is re-inserted", {
          
          set.seed(52787)
          require(BDAepimodel)
          
          popsize <- 3
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
          }
          
          
          epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.5),
                                    popsize = popsize,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = rnorm(1, 1.5, 1e-6), mu = rnorm(1, 1, 1e-6), rho = 0.5,  S0 = 0.5, I0 = 0.5, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          epimodel <- prepare_epimodel(epimodel)
          
          # array of tpms
          tpmSeqs(tpms = epimodel$tpms, pop_mat = epimodel$pop_mat, eigen_vals = epimodel$eigen_values, eigen_vecs = epimodel$eigen_vectors, inverse_vecs = epimodel$inv_eigen_vectors, irm_keys = epimodel$keys)
          
          # TPM product subsequences
          tpmProdSeqs(tpm_prods = epimodel$tpm_products, tpms = epimodel$tpms, obs_time_inds = epimodel$obs_time_inds)
          
          # Emission matrix
          epimodel$emission_mat <- build_emission_mat(emission_mat = epimodel$emission_mat, epimodel = epimodel)
          
          # FB matrices
          buildFBMats(epimodel$fb_mats, epimodel$tpm_products, epimodel$emission_mat, epimodel$initdist, epimodel$obs_time_inds)      
          
          .pop_mat <- unname(epimodel$pop_mat[1:10,])
          .config_mat <- unname(epimodel$config_mat[1:10,])
          
          # remove trajectory from the counts in config_mat 
          # and obs_mat, and update the tpm sequences to 
          # reflect the removal. 
          epimodel <- remove_trajectory(epimodel, subject = 1, save_path = TRUE)
          
          # re-insert the trajectory
          epimodel <- reinsert_path(epimodel, subject = 1)

          expect_equal(unname(epimodel$pop_mat[complete.cases(epimodel$pop_mat),]), .pop_mat[complete.cases(.pop_mat),])
})