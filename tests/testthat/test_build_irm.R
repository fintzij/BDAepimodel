library(BDAepimodel)

context("Constructing the the rate matrices")


test_that("The rate matrices are all computed correctly", {
          
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
                                    params = c(beta = 0.9, mu = 1, rho = 0.5,  S0 = 0.5, I0 = 0.5, R0 = 0), 
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
          
          expect_equal(.epimodel$.irm[["1"]][2,2], -1)
          expect_equal(.epimodel$.irm[["2"]][1,2], 1.8)
          expect_equal(.epimodel$.irm[["0"]][1,], rep(0, 3))
          expect_equal(.epimodel$.irm[["3"]][3,], rep(0, 3))
          
          # update the IRMs and retest
          build_irm(.epimodel)
          
          expect_equal(.epimodel$.irm[["1"]][1,1], -0.9)
          expect_equal(.epimodel$.irm[["2"]][2,3], 1)
          expect_equal(.epimodel$.irm[["0"]][1,], rep(0, 3))
          expect_equal(.epimodel$.irm[["3"]][3,], rep(0, 3))
}) 