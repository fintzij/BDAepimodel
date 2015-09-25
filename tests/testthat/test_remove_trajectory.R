library(BDAepimodel) 

context("Removing the contribution of a trajectory to the compartment counts")

test_that("The updated counts are correct in the configuration matrix", {
          
          set.seed(52787)
          
          popsize <- 5
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
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
          
          remove_trajectory(.epimodel, ".X2")
          
          # subject 2 is susceptible in the interval [0, 1.8705), infected in
          # the interval [1.8705, 4.44065), and recovered after. Therefore, the
          # count of susceptibles at times 0, 0.5, 1, and 1.5 should be reduced by
          # one, the count of infectives at times 2, 2.5, 3, 3.5, and 4 should
          # be reduced by one, and the counts of recovered individuals at times
          # 4.5 and 5 should be reduced by one.
          
          expect_equal(.epimodel$config_mat[which(.epimodel$config_mat[,"time"] < .epimodel$config_mat[which(.epimodel$config_mat[,"ID"] == 2)[1],"time"]) ,"S"], c(1, 1, 0, 0, 0, 0, 0, 0))
})

test_that("The observation matrix is updated correctly when a trajectory is removed", {
          
          set.seed(52787)
          
          popsize <- 5
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
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
          
          remove_trajectory(.epimodel, ".X2")
          
          # subject 2 is susceptible in the interval [0, 1.8705), infected in
          # the interval [1.8705, 4.44065), and recovered after. Therefore, the
          # count of susceptibles at times 0, 0.5, 1, and 1.5 should be reduced by
          # one, the count of infectives at times 2, 2.5, 3, 3.5, and 4 should
          # be reduced by one, and the counts of recovered individuals at times
          # 4.5 and 5 should be reduced by one.
          
          expect_equal(.epimodel$obs_mat[,"I_augmented"], c(3, 2, 1, 1, 1, 1, 0, 0, 0, 0))
})