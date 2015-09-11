# Batch file to execute test_simulate_epimodel

library(batch)
k=0
for (sim_num in 1:2) {
          k = k+1
          rbatch("./test_simulate_epimodel.R", seed = k, sim_num = sim_num)
}

