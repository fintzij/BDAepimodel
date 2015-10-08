library(BDAepimodel)

context("Checking and instatiating a new rate matrix and eigen decomposition")

test_that("It is detected whether a missing rate matrix is needed", {
          
          set.seed(52787)
          
          popsize <- 5
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
          }
          
          
          # R0 = 4, mu = 1, rho = 0.5, p0 = 0.05
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
          
          # remove the rate matrix and decompositions for keys 0 and 3
          .epimodel$.irm_keys <- .epimodel$.irm_keys[-charmatch(c(0, 3), .epimodel$.irm_keys)] 
          rm(list = c("0", "3"), envir = .epimodel$.irm)
          rm(list = c("0", "3"), envir = .epimodel$.eigen)
          
          .unique_configs <- unique(.epimodel$config_mat[1:.epimodel$.ind_final_config, .epimodel$index_states, drop = FALSE])
          .keys <- apply(.unique_configs, 1, paste0, collapse = .epimodel$index_states)
          
          # indices of .keys vector for which keys are missing
          .missing_inds <- which(is.na(charmatch(.keys, .epimodel$.irm_keys)))
          .missing_keys <- .keys[.missing_inds]
          
          expect_equal(.missing_keys, c("3", "0"))
})

test_that("Missing rate matrices and eigen decompositions are instatiated", {
          
          set.seed(52787)
          
          popsize <- 5
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          d_meas_process <- function(state, meas_vars, params, log = TRUE) {
                    dbinom(x = state[, paste(meas_vars, "_observed", sep="")], size = state[, paste(meas_vars, "_augmented", sep = "")], prob = params["rho"], log = log) 
          }
          
                    # R0 = 4, mu = 1, rho = 0.5, p0 = 0.05
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
          
          # remove the rate matrix and decompositions for keys 0 and 3
          .epimodel$.irm_keys <- .epimodel$.irm_keys[-charmatch(c(0, 3), .epimodel$.irm_keys)] 
          .epimodel$.irm_key_lookup <- .epimodel$.irm_key_lookup[- which(rownames(.epimodel$.irm_key_lookup) %in% c("0","3")), , drop = FALSE] 
          rm(list = c("0", "3"), envir = .epimodel$.irm) 
          rm(list = c("0", "3"), envir = .epimodel$.eigen) 
          
          # run check_irm to re-instatiate the required objects
          check_irm(.epimodel)
          
          .unique_configs <- unique(.epimodel$config_mat[1:.epimodel$.ind_final_config, .epimodel$index_states, drop = FALSE])
          .keys <- apply(.unique_configs, 1, paste0, collapse = .epimodel$index_states)
          
          # indices of .keys vector for which keys are missing
          .missing_inds <- which(is.na(charmatch(.keys, .epimodel$.irm_keys)))

          expect_equal(length(.missing_inds), 0)
})