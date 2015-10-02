library(BDAepimodel)

context("Calculation of the population-level log-likelihood and measurement processes log-likelihood") 

test_that("The log-likelihoods are computed correctly", {
          
          set.seed(52787)
          
          popsize <- 3
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
          }
          
          # evaluates initial distribution for a single subject
          d_initdist <- function(state, params, log = TRUE) {
                    if(log == TRUE) {
                              (state == 2)* log(params["p0"]) + (state == 1) * log(1-params["p0"]) + ifelse(state == 3, -Inf, 0)
                    } else {
                              (params["p0"] ^ (state == 2)) * ((1-params["p0"])^(state == 1)) * (0^(state == 3))
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
                                    params = c(beta = rnorm(1, 1, 1e-6), mu = rnorm(1,1,1e-6), rho = 0.5, p0 = 0.5), 
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
          
          build_irm(.epimodel)
          
          config_mat <- .epimodel$config_mat[c(1, which(.epimodel$config_mat[,"ID"] != 0), .epimodel$.ind_final_config),]
          
          pop_likelihood <- 2 * log(.epimodel$params["p0"]) + log( 1- .epimodel$params["p0"]) + 
                    log(2 * epimodel$params["beta"]) - (2 + 2)*(config_mat[2, "time"] - config_mat[1, "time"]) +
                    log(1) - (3)*(config_mat[3, "time"] - config_mat[2, "time"]) +
                    log(1) - (2)*(config_mat[4, "time"] - config_mat[3, "time"]) +
                    log(1) - (1)*(config_mat[5, "time"] - config_mat[4, "time"]) 
          
          obs_likelihood <-sum(dbinom(.epimodel$obs_mat[,"I_observed"], .epimodel$obs_mat[,"I_augmented"], prob = .epimodel$params["rho"], log= TRUE))
          
          subj1_likelihood <- log(1 - .epimodel$params["p0"]) + log(2*.epimodel$params["beta"]) - 2*.epimodel$params["beta"] * (config_mat[2,"time"] - config_mat[1,"time"]) + log(.epimodel$params["mu"])  - .epimodel$params["mu"] *(config_mat[3,"time"] - config_mat[2,"time"])
          
          expect_equal(signif(calc_pop_likelihood(epimodel = .epimodel, log = TRUE), digits = 6), as.numeric(signif(pop_likelihood, digits = 6)))
          expect_equal(as.numeric(calc_obs_likelihood(epimodel = .epimodel, log = TRUE)), as.numeric(obs_likelihood))
          expect_equal(as.numeric(calc_subj_likelihood(epimodel = .epimodel, subject = 1, subj_ID = ".X1", log = TRUE)), as.numeric(subj1_likelihood))
}) 