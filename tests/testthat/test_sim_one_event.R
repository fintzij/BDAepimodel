library(BDAepimodel)

context("Simulating the next event within simulate_epimodel")

# Tests for extended state space
test_that("Rates are computed properly", {
          # initialize epimodel
          epimodel <- init_epimodel(obstimes = 0:10,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = 0.5, mu = 1, rho = 0.5, S0 = 0.95, I0 = 0.05, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T))
          
          init_state = c(S = 7, I = 3, R = 0)
          epimodel$popsize <- sum(init_state)
          
          epimodel$config_mat <- init_config_mat(init_state = init_state, t0 = min(epimodel$obstimes), tmax = max(epimodel$obstimes))
          
          # set up items for config_list
          config <- epimodel$config_mat[epimodel$config_mat[,"time"] == min(epimodel$obstimes), , drop = FALSE]
          rate_mat <- matrix(1, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize rate matrix
          dt_mat <- matrix(Inf, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize dt matrix 
          
          config_list <- list(config = config, rate_mat = rate_mat, dt_mat = dt_mat, keep_going = TRUE)
          
          # simulate one step
          set.seed(52787)
          config_list <- sim_one_event(config_list, epimodel, lump = FALSE)
          
          expect_equal(config_list$rate_mat[1:7, 1], rep(0.5 * 3, 7)) # infectivity rates for susceptibles
          expect_equal(config_list$rate_mat[8:10, 1], rep(0, 3)) # infectivity rates for current susceptibles
          expect_equal(config_list$rate_mat[1:7, 2], rep(0, 7)) # recovery rate for susceptibles
          expect_equal(config_list$rate_mat[8:10, 2], rep(1, 3)) # recovery rate for infecteds
})

test_that("The next configuration is recorded correctly", {
          epimodel <- init_epimodel(obstimes = 0:10,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = 0.5, mu = 1, rho = 0.5, S0 = 0.95, I0 = 0.05, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T))
          
          init_state = c(S = 7, I = 3, R = 0)
          epimodel$popsize <- sum(init_state)
          
          epimodel$config_mat <- init_config_mat(init_state = init_state, t0 = min(epimodel$obstimes), tmax = max(epimodel$obstimes))
          
          # set up items for config_list
          config <- epimodel$config_mat[epimodel$config_mat[,"time"] == min(epimodel$obstimes), , drop = FALSE]
          rate_mat <- matrix(1, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize rate matrix
          dt_mat <- matrix(Inf, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize dt matrix 
          
          config_list <- list(config = config, rate_mat = rate_mat, dt_mat = dt_mat, keep_going = TRUE)
          
          # simulate one step
          set.seed(52787)
          config_list <- sim_one_event(config_list, epimodel, lump = FALSE)
          
          # the next reaction should be an infection by subject 5 at time 0.03046
          expect_equal(as.numeric(config_list$config[,"ID"]), 5) # the event is associated with subject 5
          expect_equal(as.numeric(config_list$config[,"Event"]), 1) # the event is an infection event
          expect_equal(as.numeric(config_list$config[,"time"]), min(config_list$dt_mat))
          
          # the next event is a recovery of subject 5. Check this too.
          t <- as.numeric(config_list$config[,"time"])
          
          config_list <- sim_one_event(config_list, epimodel, lump = FALSE)
          
          # the next reaction should be an infection by subject 5 at time 0.03046
          expect_equal(as.numeric(config_list$config[,"ID"]), 5) # the event is associated with subject 5
          expect_equal(as.numeric(config_list$config[,"Event"]), 2) # the event is an infection event
          expect_equal(as.numeric(config_list$config[,"time"]), t+ min(config_list$dt_mat))
})

test_that("The compartment counts are lumped correctly", {
          epimodel <- init_epimodel(obstimes = 0:10,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = 0.5, mu = 1, rho = 0.5, S0 = 0.95, I0 = 0.05, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T))
          
          init_state = c(S = 7, I = 3, R = 0)
          epimodel$popsize <- sum(init_state)
          
          epimodel$config_mat <- init_config_mat(init_state = init_state, t0 = min(epimodel$obstimes), tmax = max(epimodel$obstimes))
          
          
          # set up items for config_list
          config <- epimodel$config_mat[epimodel$config_mat[,"time"] == min(epimodel$obstimes), , drop = FALSE]
          rate_mat <- matrix(1, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize rate matrix
          dt_mat <- matrix(Inf, nrow = epimodel$popsize, ncol = nrow(epimodel$flow)) # initialize dt matrix 
          
          config_list <- list(config = config, rate_mat = rate_mat, dt_mat = dt_mat, keep_going = TRUE)
          
          # simulate one step
          set.seed(52787)
          config_list <- sim_one_event(config_list, epimodel, lump = FALSE)
          
          # the next reaction should be an infection by subject 5 at time 0.03046
          expect_equal(sum(config_list$config[,7:16] == 1), 6) # there should be six susceptibles
          expect_equal(sum(config_list$config[,7:16] == 2), 4) # 4 infecteds
          expect_equal(sum(config_list$config[,7:16] == 3), 0) # no recovereds
          
          # the next infection is a recovery, check that too
          config_list <- sim_one_event(config_list, epimodel, lump = FALSE)
          expect_equal(sum(config_list$config[,7:16] == 1), 6) # there should be six susceptibles
          expect_equal(sum(config_list$config[,7:16] == 2), 3) # 4 infecteds
          expect_equal(sum(config_list$config[,7:16] == 3), 1) # no recovereds
})

# Test for the lumped state space
test_that("The lumped rates are computed properly", {
          
          # initialize the epimodel object
          epimodel <- init_epimodel(obstimes = 0:10,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = 0.5, mu = 1, rho = 0.5, S0 = 0.95, I0 = 0.05, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T))
          
          init_state <- c(S = 7, I = 3, R = 0)
          epimodel$popsize <- sum(init_state)
          epimodel$config_mat <- init_config_mat(init_state = init_state, t0 = min(epimodel$obstimes), tmax = max(epimodel$obstimes))
          
          #set up the items for config_list
          config <- epimodel$config_mat[epimodel$config_mat[,"time"] == min(epimodel$obstimes), , drop = FALSE]
          rate_mat <- matrix(1, nrow = 1, ncol = nrow(epimodel$flow)) # initialize rate matrix
          dt_mat <- matrix(Inf, nrow = 1, ncol = nrow(epimodel$flow)) # initialize dt matrix
          
          config_list <- list(config = config, rate_mat = rate_mat, dt_mat = dt_mat, keep_going = TRUE)
          
          # simulate one step
          set.seed(52787)
          config_list <- sim_one_event(config_list, epimodel, lump = TRUE)
          
          expect_equal(as.numeric(config_list$rate_mat), c(0.5 * 7 * 3, 1 * 3))
})

test_that("Subjects are always chosen from the appropriate risk set", {
          
          # initialize the epimodel object
          epimodel <- init_epimodel(obstimes = 0:10,
                                    states = c("S", "I", "R"), 
                                    params = c(beta = 0.5, mu = 1, rho = 0.5, S0 = 0.95, I0 = 0.05, R0 = 0), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T))
          
          init_state <- c(S = 7, I = 3, R = 0)
          epimodel$popsize <- sum(init_state)
          epimodel$config_mat <- init_config_mat(init_state = init_state, t0 = min(epimodel$obstimes), tmax = max(epimodel$obstimes))
          
          #set up the items for config_list
          config <- epimodel$config_mat[epimodel$config_mat[,"time"] == min(epimodel$obstimes), , drop = FALSE]
          rate_mat <- matrix(1, nrow = 1, ncol = nrow(epimodel$flow)) # initialize rate matrix
          dt_mat <- matrix(Inf, nrow = 1, ncol = nrow(epimodel$flow)) # initialize dt matrix
          
          config_list <- list(config = config, rate_mat = rate_mat, dt_mat = dt_mat, keep_going = TRUE)
          
          # redraw the next step many times
          set.seed(52787)
          next_events <- replicate(10000, sim_one_event(config_list, epimodel, lump = TRUE)$config[,c("ID", "Event")])
          ID_group <- ifelse(next_events["ID",] <8, 1, 2)
          
          expect_equal(ID_group, next_events["Event", ])
})