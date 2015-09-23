library(BDAepimodel)

context("Calculation of the population-level log-likelihood and measurement processes log-likelihood") 

test_that("The population-level log-likelihood is computed correctly", {
          
          set.seed(52787)
          
          popsize <- 3
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_truth", sep = "")], prob = params["rho"], log = log) 
          }
          
          # evaluates initial distribution for a single subject
          d_initdist <- function(state, params, log = TRUE) {
                    if(log == TRUE) {
                              (state == 2)* log(params["p0"]) + (state == 1) * log(1-params["p0"])
                    } else {
                              (params["p0"] ^ (state == 2)) * ((1-params["p0"])^(state == 1))
                    }
          }
          
          # subject level simulation of initial state at time t0
          r_initdist <- function(params) {
                    sample.int(3, 1, prob = c(1-params["p0"], params["p0"], 0))
          }
          
          # R0 = 4, mu = 1, rho = 0.5, p0 = 0.05
          epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.5),
                                    popsize = popsize,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = 1, mu = 1, rho = 0.5, p0 = 0.5), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process,
                                    d_initdist = d_initdist,
                                    r_initdist = r_initdist)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          .epimodel <- list2env(epimodel, parent = emptyenv(), hash = TRUE)
          
          .epimodel$.config_inds <- which(grepl(".X", colnames(.epimodel$config_mat)))
          
          expand_config_mat(.epimodel)
          
          config_mat <- .epimodel$config_mat[c(1, which(.epimodel$config_mat[,"ID"] != 0), .epimodel$.ind_final_config),]
          
          pop_likelihood <- log(.epimodel$params["p0"]) + 2 * log( 1- .epimodel$params["p0"]) + 
                    log(1) - (1 + 2)*(config_mat[2, "time"] - config_mat[1, "time"]) +
                    log(2) - (2 + 2)*(config_mat[3, "time"] - config_mat[2, "time"]) +
                    log(1) - (3)*(config_mat[4, "time"] - config_mat[3, "time"]) +
                    log(1) - (2)*(config_mat[5, "time"] - config_mat[4, "time"]) +
                    log(1) - (1)*(config_mat[6, "time"] - config_mat[5, "time"]) -
                    (0) * (config_mat[7, "time"] - .epimodel$config_mat[6, "time"]) 
          
          obs_likelihood <-sum(dbinom(.epimodel$obs_mat[,"I_observed"], .epimodel$obs_mat[,"I_truth"], prob = .epimodel$params["rho"], log= TRUE))
          
          expect_equal(as.numeric(calc_likelihoods(epimodel = .epimodel)[[1]]), as.numeric(pop_likelihood))
          expect_equal(as.numeric(calc_likelihoods(epimodel = .epimodel)[[2]]), as.numeric(obs_likelihood))
}) 