# Simulation test to check that simulate_epimodel produces the same results as 
# the GillespieSSA package. Two models are simulated: SIR and SIRS.
# Each of the simulations will produce a plot comparing the trajectories
# produced by simulate epimodel with the trajectories produced by GillespieSSA.

library(BDAepimodel)
library(GillespieSSA)
library(doParallel)
library(itertools)
library(foreach)
library(iterators)
library(ggplot2)
library(reshape2)
library(batch)

parseCommandArgs()

# Set up simulation objects ----------------------------------------------------------------

niter <- 50; nobs <- length(seq(0, 10, by = 0.05))

if(sim_num == 1){
          # Set up BDAepimodel objects
          initialization_fcn <- function(){
                    init_state <- as.numeric(rmultinom(1, 200, c(0.95, 0.05, 0)))
                    names(init_state ) <- c("S", "I", "R")
                    return(init_state)
          }
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          # R0 = 4, mu = 1, rho = 0.5, p0 = 0.05
          epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.05),
                                    states = c("S", "I", "R"), 
                                    params = c(beta = 0.02, mu = 1, rho = 0.5, p0 = 0.05), 
                                    rates = c("beta * I", "mu"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process)
          
          # Set up GillespieSSA objects
          params <- epimodel$params
          a = c("beta*S*I", "mu*I")
          nu <- matrix(c(-1,+1,0, 0, -1, +1),nrow=3)
          gSSA_traj <- rep(0, nobs)
          
} else if(sim_num == 2) {
          results_SIRS <- data.frame(method = rep(c("BDAepimodel", "GillespieSSA"), each = nobs),
                                    estimate = rep(0, 2*nobs),
                                    lower = rep(0, 2*nobs),
                                    upper = rep(0, 2*nobs))
          
          # Set up BDAepimodel objects
          initialization_fcn <- function(){
                    init_state <- as.numeric(rmultinom(1, 200, c(0.95, 0.05, 0)))
                    names(init_state ) <- c("S", "I", "R")
                    return(init_state)
          }
          
          r_meas_process <- function(state, meas_vars, params){
                    rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          }
          
          # R0 = 4, mu = 1, rho = 0.5, p0 = 0.05
          epimodel <- init_epimodel(obstimes = seq(0, 10, by = 0.05),
                                    states = c("S", "I", "R"), 
                                    params = c(beta = 0.02, mu = 1, lambda = 0.5, rho = 0.5, p0 = 0.05), 
                                    rates = c("beta * I", "mu", "lambda"), 
                                    flow = matrix(c(-1, 1, 0, 0, -1, 1, 1, 0, -1), ncol = 3, byrow = T), 
                                    meas_vars = "I",
                                    r_meas_process = r_meas_process)
          
          # Set up GillespieSSA objects
          params <- epimodel$params
          a = c("beta*S*I", "mu*I", "lambda * R")
          nu <- matrix(c(-1,+1,0, 0, -1, +1, 1, 0, -1),nrow=3)
          gSSA_traj <- rep(0, nobs)
}


# Register the cluster and run the simulation -----------------------------

cl <- makeForkCluster(15)
cl <- makeCluster(3)
registerDoMC(cl)

trajecs <- foreach(j = 1:niter, .combine = 'cbind', .packages = c("BDAepimodel", "GillespieSSA")) %dopar% {
          
          init_state <- initialization_fcn()
          if(init_state["I"] == 0){
                    while(init_state["I"] == 0){
                              init_state <- initialization_fcn()
                    }
          }
          
          BDAepimodel_traj <- simulate_epimodel(epimodel = epimodel, init_state = init_state)$obs_mat[,"I_truth"] 
          
          gSSA_sim <- ssa(init_state, a, nu, params, tf = 10)$data[,c(1,3)]
          
          for(k in 1:length(gSSA_traj)) {
                    gSSA_traj[k] <- gSSA_sim[sum(gSSA_sim[,1] <= epimodel$obstimes[k]), 2]
          }
          
          return(c(BDAepimodel_traj, gSSA_traj))
} 


# Process results and generate plots --------------------------------------

results <- data.frame(method = rep(c("BDAepimodel", "GillespieSSA"), each = nobs),
                      time = seq(0, 10, by = 0.05),
                      estimate = rep(0, 2*nobs),
                      lower = rep(0, 2*nobs),
                      upper = rep(0, 2*nobs))

results$estimate <- rowMeans(trajecs)
results$lower <- pmax(0,results$estimate - 1.96 * apply(trajecs, 1, sd) / sqrt(niter))
results$upper <- results$estimate + 1.96 * apply(trajecs, 1, sd) / sqrt(niter)

pdf("test_simulate_epimodel.pdf")

ggplot(results, aes(x = time, y = estimate, colour = method)) + geom_line() + geom_ribbon(data = results, aes(x = time, ymin = lower, ymax = upper, fill = method), alpha = 0.1) 

dev.off()