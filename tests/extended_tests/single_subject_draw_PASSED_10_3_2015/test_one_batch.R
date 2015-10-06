# Batch file to execute test_simulate_epimodel

library(batch)
k=0
for (sim_num in 1) {
          k = k+1
          rbatch("./BDAepimodel_draw_one.R", seed = k)
}

