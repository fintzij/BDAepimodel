# Simulation test to check that simulate_epimodel produces the same results as 
# the GillespieSSA package. Two models are simulated: SIR and SIRS.
# Each of the simulations will produce a plot comparing the trajectories
# produced by simulate epimodel with the trajectories produced by GillespieSSA.

library(BDAepimodel)
library(GillespieSSA)
library(doParallel)
library(doMC)
library(itertools)
library(foreach)
library(iterators)
library(ggplot2)
library(reshape2)
library(batch)

parseCommandArgs()

# Set up simulation objects ----------------------------------------------------------------
# Set up GillespieSSA objects
niter <- 10000
popsize = 10
obstimes = seq(0, 10, by = 1)
nobs <- length(obstimes)
params <- c(beta = 0.5, mu = 1, rho = 0.5, S0 = 0.8, I0 = 0.2, R0 = 0)
a = c("beta*S*I", "mu*I")
nu <- matrix(c(-1,+1,0, 0, -1, +1),nrow=3)

initialization_fcn <- function(){
          init_state <- as.numeric(rmultinom(1, popsize, c(0.8, 0.2, 0)))
          names(init_state ) <- c("S", "I", "R")
          return(init_state)
}

r_meas_process <- function(state, meas_vars, params){
          rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
}

epimodel <- init_epimodel(obstimes = obstimes,
                          popsize = popsize,
                          states = c("S", "I", "R"), 
                          params = params, 
                          rates = c("beta * I", "mu"), 
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                          meas_vars = "I",
                          r_meas_process = r_meas_process)


# Register the cluster and run the simulation -----------------------------

# nclust = 15
# cl <- makeForkCluster(nclust)

trajecs <- matrix(0, nrow = 3*length(obstimes), ncol = niter)
simSIR_step1 <- rep(0, niter)
gssa_step1 <- rep(0, niter)

for(j in 1:niter) {
          
          if(j%%100 == 0) print(j)
          
          init_state <- initialization_fcn()
          if(init_state["I"] == 0){
                    while(init_state["I"] == 0){
                              init_state <- initialization_fcn()
                    }
          }
          
          simulateSIR_traj <- simulate_SIR(obstimes = obstimes, params = params, popsize = popsize, init_config = init_state, trim = FALSE)$obs_mat[,"I_augmented"]

          BDAepimodel_traj <- simulate_epimodel(epimodel = epimodel, init_state = init_state, lump = TRUE, trim = FALSE)$obs_mat[,"I_augmented"]
          
          gSSA_sim <- ssa(init_state, a, nu, params, tf = 10)$data[,c(1,3)]
          
          gSSA_traj <- rep(0, nobs)
          for(k in 1:length(gSSA_traj)) {
                    gSSA_traj[k] <- gSSA_sim[sum(gSSA_sim[,1] <= obstimes[k]), 2]
          }
          
          trajecs[,j] <- c(simulateSIR_traj, BDAepimodel_traj, gSSA_traj)
}

# Process results and generate plots --------------------------------------

results <- data.frame(method = rep(c("simulateSIR","BDAepimodel", "GillespieSSA"), each = nobs),
                      time = obstimes,
                      estimate = rep(0, 3*nobs),
                      lower = rep(0, 3*nobs),
                      upper = rep(0, 3*nobs))

results$estimate <- rowMeans(trajecs)
results$lower <- pmax(0,results$estimate - 1.96 * apply(trajecs, 1, sd) / sqrt(niter))
results$upper <- results$estimate + 1.96 * apply(trajecs, 1, sd) / sqrt(niter)

results_diff <- data.frame(time = obstimes,
                           estimate = rep(0, nobs),
                           lower = rep(0, nobs), 
                           upper = rep(0, nobs))
results_diff$estimate <- rowMeans(trajecs[1:nobs, ]) - rowMeans(trajecs[(nobs+1):(2*nobs), ])
results_diff$lower <- results_diff$estimate - 1.96 * (apply(trajecs[1:nobs, ], 1, sd) + apply(trajecs[(nobs+1):(2*nobs), ], 1, sd)) / sqrt(niter)
results_diff$upper <- results_diff$estimate + 1.96 * (apply(trajecs[1:nobs, ], 1, sd) + apply(trajecs[(nobs+1):(2*nobs), ], 1, sd)) / sqrt(niter)

pdf(paste("test_simulate_epimodel_extendedSS",sim_num,".pdf", sep = ""))

ggplot(results, aes(x = time, y = estimate, colour = method)) + geom_line(size = 1) + geom_ribbon(data = results, aes(x = time, ymin = lower, ymax = upper, fill = method), alpha = 0.1, size = 0) + labs(y = "Average number of infecteds") + theme_bw()

ggplot(results_diff, aes(x = time, y = estimate)) + geom_line(size = 1) + geom_ribbon(data = results_diff, aes(x = time, ymin = lower, ymax = upper), alpha = 0.1, size = 0) + labs(y = "Average difference in number of infecteds") + theme_bw()

dev.off()