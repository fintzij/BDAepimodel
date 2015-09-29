library(BDAepimodel)

context("Construction of the forward-backward matrices")

test_that("FB matrix construction in BDAepimodel matches the tested augSIR build", {
          
          # instatiate the augSIR function
          fwd_augSIR <- function(tpms, emits, initdist){
                    fbmats <- array(0, dim = dim(tpms))
                    
                    pi0 <- normalize(initdist * emits[,1])
                    fbmats[,,1] <- normalize(outer(pi0,emits[,2], FUN="*") * tpms[,,1])
                    
                    for(s in 2:dim(fbmats)[3]){
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
                                    params = c(beta = rnorm(1, 0.5, 1e-6), mu = rnorm(1, 1, 1e-6), rho = 0.5, p0 = 0.2), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process,
                                    d_meas_process = d_meas_process,
                                    d_initdist = d_initdist,
                                    r_initdist = r_initdist)
          
          epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = FALSE)
          
          # convert epimodel list to environment to enable multiple assignment
          .epimodel <- list2env(epimodel, parent = emptyenv(), hash = TRUE)
          
          # identify the unmeasured compartments, constants for number of
          # measured and unmeasured states
          .epimodel$unmeasured_vars     <- setdiff(.epimodel$states, .epimodel$meas_vars)
          .epimodel$num_unmeasured      <- length(.epimodel$unmeasured_vars)
          .epimodel$num_measured        <- length(.epimodel$meas_vars)
          .epimodel$num_states          <- length(.epimodel$states)
          
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
          
          # reset the vector of tpms to build (set all to false at the end of the build), then build the tpm
          # sequences
          .epimodel$.tpms_to_build <- get_tpms_to_build(.epimodel)
          build_tpm_seqs(.epimodel)
          
          # build the emission matrix
          build_emission_mat(.epimodel)
          
          # construct the forward-backward matrices
          build_fb_mats(.epimodel)
          
          augSIR_tpms <- array(0, dim = c(3,3,.epimodel$nobs - 1))
          for(k in 1:dim(augSIR_tpms)[3]){
                    augSIR_tpms[,,k] <- .epimodel$.tpm_products[[which(.epimodel$config_mat[,"ID"] == 0)[k]]]
          }
          
          augSIR_fbarray <- fwd_augSIR(tpms = augSIR_tpms, emits = .epimodel$.emission_mat, initdist = .epimodel$.initdist)
          
          expect_true(all(simplify2array(.epimodel$.fb_mats) == augSIR_fbarray))
})