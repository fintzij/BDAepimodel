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
                              (state == 2)* log(params["p0"]) + (state == 1) * log(1-params["p0"]) + (state == 3) * log(0)
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

test_that("The sequence of tpm products is computed correctly when a trajectory is removed", {
          
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
                                    params = c(beta = rnorm(1, 1, 1e-6), mu = rnorm(1, 1, 1e-6), rho = 0.5, p0 = 0.5), 
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
          
          .epimodel$.tpms_to_build <- get_tpms_to_build(.epimodel)
          
          build_tpm_seqs(.epimodel)
          
          remove_trajectory(.epimodel, 1)
          
          # check to see if any additional irms are needed.
          # if so, check_irm will instatiate the required
          # matrices and their eigen decompositions
          check_irm(.epimodel)
          
          # get indices of tpms that need to be rebuilt
          .epimodel$.tpms_to_build <- get_tpms_to_build(.epimodel, subject = 1)
          
          build_tpm_seqs(.epimodel)
          
          .tpm1 <-  eigen(.epimodel$.irm[["4"]])$vectors %*% 
                    diag(exp(eigen(.epimodel$.irm[["4"]])$values * (.epimodel$config_mat[2,"time"] - .epimodel$config_mat[1,"time"]))) %*% solve(eigen(.epimodel$.irm[["4"]])$vectors)
          
          .tpm2 <- eigen(.epimodel$.irm[["3"]])$vectors %*% 
                    diag(exp(eigen(.epimodel$.irm[["3"]])$values * (.epimodel$config_mat[3,"time"] - .epimodel$config_mat[2,"time"]))) %*% solve(eigen(.epimodel$.irm[["3"]])$vectors)
          
          expect_equal(.epimodel$.tpms[[2]], .tpm2)
          expect_equal(.epimodel$.tpm_products[[1]], .tpm1 %*% .tpm2)
})

test_that("tpms are computed correctly following multiple removals and re-insertions", {
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
                                    params = c(beta = rnorm(1, 1, 1e-6), mu = rnorm(1, 1, 1e-6), rho = 0.5, p0 = 0.5), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process,
                                    d_initdist = d_initdist,
                                    r_initdist = r_initdist)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
          
          .epimodel <- list2env(epimodel, parent = emptyenv(), hash = TRUE)
          
          
          # initialize the list of log-likelihoods
          epimodel$likelihoods <- list(pop_likelihood_cur = NULL,
                                       pop_likelihood_new = NULL,
                                       subj_likelihood_cur = NULL,
                                       subj_likelihood_new = NULL,
                                       obs_likelihood  = NULL)
          
          # identify the unmeasured compartments, constants for number of
          # measured and unmeasured states
          .epimodel$unmeasured_vars     <- setdiff(.epimodel$states, .epimodel$meas_vars)
          .epimodel$num_unmeasured      <- length(.epimodel$unmeasured_vars)
          .epimodel$num_measured        <- length(.epimodel$meas_vars)
          .epimodel$num_states          <- length(.epimodel$states)
          
          # identify whether there is structure in the flow matrix that can be
          # leveraged when resampling subject-level trajectories
          detect_structure(.epimodel)
          

          # add buffer to the configuration matrix, also adds configurations at
          # observation times and instatiates .epimodel$.ind_final_config
          .epimodel$config_mat <- expand_config_mat(.epimodel)
          
          # Initialize two lists for the transition probability matrices, one 
          # for the matrices and one for products. Also initialize a bookkeeping
          # vector indicating which tpms and tpm product matrices need to be 
          # recomputed. All tpms and products need to be recomputed following 
          # parameter updates, but only the tpms and products for the intervals 
          # affected by a change in one of the compartments indicated by 
          # index_states need to be recomputed when a new trajectory is redrawn.
          .epimodel$.tpms               <- vector("list", length = nrow(.epimodel$config_mat))
          .epimodel$.tpm_products       <- vector("list", length = nrow(.epimodel$config_mat))
          
          # initialize the matrix of emission probabilities
          .epimodel$nobs                <- length(.epimodel$obstimes)
          .epimodel$.emission_mat       <- matrix(0, nrow = .epimodel$num_states, ncol = .epimodel$nobs, dimnames = list(.epimodel$states, .epimodel$obstimes))
          
          # initialize the forward-backward matrices
          .epimodel$.fb_mats             <- vector("list", length = .epimodel$nobs - 1)
          
          # initialize the vector of subject-level initial state probabilities
          .epimodel$.initdist <- build_initdist(.epimodel)
          
          # get indices for the subject configuration portion of the configuration matrix
          .epimodel$.config_inds <- which(grepl(".X", colnames(.epimodel$config_mat)))
          
          build_irm(.epimodel)
          
          for(k in 1:3) {
                    
                    # re-compute the array of rate matrices - only needs to be
                    # recomputed every time parameters are updated.
                    
                    # reset the vector of tpms to build (set all to false at the end of the build), then build the tpm
                    # sequences
                    .epimodel$.tpms_to_build <- get_tpms_to_build(.epimodel)
                    build_tpm_seqs(.epimodel)
                    
                    # re-compute the population level likelihood, measurement
                    # process likelihood, and complete data log-likelihood
                    .epimodel$likelihoods$pop_likelihood_cur <- calc_pop_likelihood(epimodel = .epimodel, log = TRUE)
                    .epimodel$likelihoods$obs_likelihood_cur <- calc_obs_likelihood(epimodel = .epimodel, log = TRUE)
                    
                    # choose which subjects should be redrawn
                    .subjects <- sample.int(n = .epimodel$popsize, size = 5)
                    
                    # cycle through subject-level trajectories to be re-drawn
                    for(j in 1:5) {
                              
                              # check to see if any additional irms are needed.
                              # if so, check_irm will instatiate the required
                              # matrices and their eigen decompositions
                              check_irm(.epimodel)
                              
                              # remove trajectory from the counts in config_mat 
                              # and obs_mat, and update the tpm sequences to 
                              # reflect the removal. 
                              remove_trajectory(.epimodel, subject = .subjects[j])
                              
                              # update instatiate missing IRMs, update the tpms 
                              # and tpm products, the emission probability mtx,
                              # and the FB matrices.
                              update_matrices(.epimodel, subject = .subjects[j])
                              
                              # draw a new subject level trajectory
                              draw_trajec(.epimodel, subject = .subjects[j])
                              
                    }
                    
          }
})