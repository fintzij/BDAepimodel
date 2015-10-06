library(BDAepimodel)
library(pryr)

context("Creating function lists")

# tests relate to extract_rate_fcns 
test_that("Rates are computed properly for a single system configuration", {
          rates <- c("beta * I", "mu")
          # state is a matrix argument, params is a vector argument
          state <- cbind(S = 10, I = 5, R = 0); states = colnames(state)
          params <- c(beta = 0.5, mu = 1, rho = 0.25, S0 = 0.5, I0 = 0.5, R0 = 0); param_names <- names(params)
          rate_fcns <- extract_rate_fcns(rates = rates, states = states, param_names = param_names)
          
          expect_equal(rate_fcns[[1]](state = state, params = params), 2.5) 
          expect_equal(rate_fcns[[2]](state = state, params = params), 1)
})

test_that("Rates are computed properly for a matrix of system configuations",{
          rates <- c("beta * I", "mu")
          state = cbind(S = seq(10,20,by=2), I = seq(10,20,by=2)/2, R = 0); states <- colnames(state)
          params <- c(beta = 0.5, mu = 1, rho = 0.25, S0 = 0.5, I0 = 0.5, R0 = 0); param_names <- names(params)
          
          rate_fcns <- extract_rate_fcns(rates = rates, states = states, param_names = param_names)
          
          expect_equal(rate_fcns[[1]](state = state, params = params), state[,"I"]*0.5)
})

test_that("Rates are computed properly when a rate depends on several compartments", {
          rates <- c("beta * S * I")
          state <- cbind(S = 10, I = 5, R = 0); states = colnames(state)
          params <- c(beta = 0.5, mu = 1, rho =0.25, S0 = 0.5, I0 = 0.5, R0 = 0); param_names <- names(params)
          rate_fcns <- extract_rate_fcns(rates = rates, states = states, param_names = param_names)
          
          expect_equal(rate_fcns[[1]](state = state, params = params), 10*5*0.5)
})