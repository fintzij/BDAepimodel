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
                                    params = c(beta = 0.9, mu = 1, rho = 0.5, p0 = 0.5), 
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
          
          .epimodel$.tpms_to_build      <- 1:(.epimodel$.ind_final_config - 1)
          build_tpm_seqs(.epimodel)
          
          tpm1 <- eigen(.epimodel$.irm[["3"]])$vectors %*% 
                    diag(exp(eigen(.epimodel$.irm[["3"]])$values * (.epimodel$config_mat[2, "time"] - .epimodel$config_mat[1, "time"]))) %*%
                    solve(eigen(.epimodel$.irm[["3"]])$vectors )
          
          tpm6 <- eigen(.epimodel$.irm[["1"]])$vectors %*% 
                    diag(exp(eigen(.epimodel$.irm[["1"]])$values * (.epimodel$config_mat[7, "time"] - .epimodel$config_mat[6, "time"]))) %*% solve(eigen(.epimodel$.irm[["1"]])$vectors )
          
          tpm7 <- eigen(.epimodel$.irm[["1"]])$vectors %*% 
                    diag(exp(eigen(.epimodel$.irm[["1"]])$values * (.epimodel$config_mat[8, "time"] - .epimodel$config_mat[7, "time"]))) %*% solve(eigen(.epimodel$.irm[["1"]])$vectors )
          
          expect_equal(.epimodel$.tpms[[1]], tpm1)
          expect_equal(.epimodel$.tpms[[6]], tpm6)
          expect_equal(.epimodel$.tpms[[7]], tpm7)
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
                                    params = c(beta = 0.9, mu = 1, rho = 0.5, p0 = 0.5), 
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
          
          .epimodel$.tpms_to_build      <- 1:(.epimodel$.ind_final_config - 1)
          build_tpm_seqs(.epimodel)
          
          tpm_prod9 <- .epimodel$.tpms[[9]] %*% .epimodel$.tpms[[10]]
          
          expect_equal(.epimodel$.tpm_products[[9]], tpm_prod9)
          
})

test_that("The tpm product sequence is changed appropriately when a single event time is changed", {
          
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
                                    params = c(beta = 0.9, mu = 1, rho = 0.5, p0 = 0.5), 
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
          
          .epimodel$.tpms_to_build <- 1:(.epimodel$.ind_final_config - 1)
          build_tpm_seqs(.epimodel)
          
          .epimodel$config_mat[10,"time"] <- 2.283
          
          .epimodel$.tpms_to_build <- 9:10
          build_tpm_seqs(.epimodel)
          
          tpm_prod9 <- .epimodel$.tpms[[9]] %*% .epimodel$.tpms[[10]]
          
          expect_equal(.epimodel$.tpm_products[[9]], tpm_prod9)
          
})