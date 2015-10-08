library(BDAepimodel)

context("Generating the correct sequence of tpm and tpm products after a trajectory is removed")

test_that("The correct indices for tpms to be rebuilt are retrieved", {
          
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
                                    params = c(beta = rnorm(1, 1, 1e-6), mu = rnorm(1, 1, 1e-6), rho = 0.5, S0 = 0.5, I0 = 0.5, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          .epimodel <- list2env(epimodel, parent = emptyenv(), hash = TRUE)
          
          .epimodel$.config_inds <- which(grepl(".X", colnames(.epimodel$config_mat)))
          
          expand_config_mat(.epimodel)
          
          build_irm(.epimodel)
          
          .epimodel$.initdist <- build_initdist(.epimodel)
          
          remove_trajectory(.epimodel, 1, save_path = TRUE)
          
          # check to see if any additional irms are needed.
          # if so, check_irm will instatiate the required
          # matrices and their eigen decompositions
          check_irm(.epimodel)
          
          # get indices of tpms that need to be rebuilt
          .epimodel$.tpms_to_build <- get_tpms_to_build(.epimodel, subject = 1)
          
          
          # subject 2 was infected in [0.4167, 3.7778). Rows 3 and 14 should be set to NA and the configuration matrix re-ordered. The tpms indices should be 2:12. 
          expect_true(all(2:12 %in%.epimodel$.tpms_to_build))
          expect_equal(as.numeric(.epimodel$config_mat[3,1]), 0.5)
          expect_equal(as.numeric(.epimodel$config_mat[13,1]), 4)
})

