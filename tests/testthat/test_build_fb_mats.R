library(BDAepimodel)

context("Construction of the forward-backward matrices")

test_that("FB matrix construction in BDAepimodel matches the tested augSIR build", {

          # instatiate the augSIR function
          fwd_augSIR <- function(tpms, emits, initdist){
                    fbmats <- array(0, dim = c(nrow(emits), nrow(emits), ncol(emits) - 1))

                    pi0 <- normalize(initdist * emits[,1])
                    fbmats[,,1] <- normalize(outer(pi0,emits[,2], FUN="*") * tpms[,,1])

                    for(s in 2:(dim(emits)[2] - 1)){
                              fbmats[,,s] <- normalize(outer(colSums(fbmats[,,s-1]), emits[,s+1], FUN="*") * tpms[,,s])

                    }

                    return(fbmats)
          }

          
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

          # array of tpms
          tpmSeqs(tpms = epimodel$tpms, pop_mat = epimodel$pop_mat, eigen_vals = epimodel$eigen_values, eigen_vecs = epimodel$eigen_vectors, inverse_vecs = epimodel$inv_eigen_vectors, irm_keys = epimodel$keys)
          
          # TPM product subsequences
          tpmProdSeqs(tpm_prods = epimodel$tpm_products, tpms = epimodel$tpms, obs_time_inds = epimodel$obs_time_inds)
          
          # Emission matrix
          epimodel$emission_mat <- build_emission_mat(emission_mat = epimodel$emission_mat, epimodel = epimodel)
          
          # FB matrices
          buildFBMats(epimodel$fb_mats, epimodel$tpm_products, epimodel$emission_mat, epimodel$initdist, epimodel$obs_time_inds)      

          augSIR_fbarray <- fwd_augSIR(tpms = epimodel$tpm_products[,,epimodel$obs_time_inds[1:(length(epimodel$obs_time_inds - 1))]], emits = epimodel$emission_mat, initdist = epimodel$initdist)

          expect_true(all(signif(epimodel$fb_mats, digits = 9) == signif(augSIR_fbarray, digits = 9)))
})