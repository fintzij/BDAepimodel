# library(BDAepimodel)
# 
# context("Initializing the subject-level bookkeeping matrix")
# 
# 
# # tests relate to extract_rate_fcns 
# test_that("The correct number of subjects are allocated to each compartment", {
#           epimodel <- init_epimodel(obstimes = 0:10,states = c("S", "I", "R"), params = c(beta = 0.5, mu = 1, rho = 0.5, p0 = 0.05), rates = c("beta * I", "mu"), flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = TRUE), init_state = c(S = 40, I = 5, R = 5)) 
#           
#           expect_equal(colSums(epimodel$subj_mat[,c("S", "I", "R")]), c(S = 40, I = 5, R = 5))
# })
