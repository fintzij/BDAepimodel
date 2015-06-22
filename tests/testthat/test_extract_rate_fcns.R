library(BDAepimodel)
library(pryr)

context("Creating function lists")

# tests relate to extract_rate_fcns 
test_that("Rates are computed properly for a single system configuration", {
          rates <- c("beta * S * I", "mu * I")
          compartments <- c("S", "I", "R")
          proc_param_names <- c("beta", "mu")
          rate_fcns <- extract_rate_fcns(rates = rates, compartments = compartments, params = proc_param_names)
          
          expect_equal(rate_fcns[[1]](10, 5, 0.5), 25)
          expect_equal(rate_fcns[[2]](10, 0.5), 5)
})

test_that("Rates are computed properly for a vector of system configuations",{
          rates <- c("beta * S * I", "mu * I")
          compartments <- c("S", "I", "R")
          proc_param_names <- c("beta", "mu")
          
          rate_fcns <- extract_rate_fcns(rates = rates, compartments = compartments, params = proc_param_names)

          S_vec <- seq(10,20,by=2)
          I_vec <- S_vec/2
          
          expect_equal(rate_fcns[[1]](S_vec, I_vec, 0.5), S_vec*I_vec*0.5)
          expect_equal(rate_fcns[[2]](I_vec, mu =1), I_vec)
})