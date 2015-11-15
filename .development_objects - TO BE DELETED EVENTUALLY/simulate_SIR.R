#' Simulate an epidemic with SIR dynamics and draw binomial samples.
#' 
#' @param obstimes vector of observation times
#' @param params vector of parameters, must be in the following order: beta, mu,
#'   rho, S0, I0, R0
#' @param popsize population size
#' @param trim logical for whether to trim the observations after the epidemic 
#'   has died off
#'   
#' @return list with a population bookkeeping matrix and an observation matrix
#' @export
#' 
#' @examples obstimes <- 0:10
#' params <- c(0.04, 1, 0.5, 0.9, 0.03, 0.07)
#' popsize <- 100
#' 
#' epidemic <- simulate_SIR(obstimes, params, popsize)
simulate_SIR <- function(obstimes, params, popsize, trim = TRUE) {
          
          # generate initial configuration
          init_config <- rmultinom(1, popsize, params[4:6])
          if(init_config[2] == 0) {
                    while(init_config[2] == 0) {
                              init_config <- rmultinom(1, popsize, params[4:6])
                    }
          }
          
          # simulate the pop_mat
          pop_mat <- simulateSIR(obstimes = obstimes, params = params[1:2], init_config = init_config)
          
          # get rid of the empty rows
          pop_mat <- pop_mat[-which(pop_mat[,1] == 0),]
          
          # add the observation times
          pop_mat <- rbind(pop_mat, cbind(obstimes, matrix(0, nrow = length(obstimes), ncol = 5)))
          
          # reorder the matrix
          row_ord <- order(pop_mat[,1])
          pop_mat <- reorderMat(pop_mat, row_ord)
          
          # attach column names
          colnames(pop_mat) <- c("time", "ID", "Event", "S", "I", "R")
          
          # set the compartment counts at observation times
          obs_time_inds <- which(pop_mat[,"ID"] == 0)
          pop_mat[1,c("S", "I", "R")] <- init_config
          
          for(j in 2:length(obs_time_inds)) {
                    pop_mat[obs_time_inds[j], c("S", "I", "R")] <- pop_mat[obs_time_inds[j] - 1, c("S", "I", "R")]
          }
          
          # initialize the observation matrix
          obs_mat <- matrix(0, nrow = length(obstimes), ncol = 3)
          colnames(obs_mat) <- c("time", "I_observed", "I_augmented")
          obs_mat[,"time"] <- obstimes
          obs_mat[,"I_augmented"] <- pop_mat[obs_time_inds, "I"]
          obs_mat[,"I_observed"] <- rbinom(length(obstimes), obs_mat[,"I_augmented"], params[3])
          
          # if trim is TRUE, then trim off the excess 0 measurements in the case
          # when the epidemic died off before the final observation time.
          if((trim == TRUE) && any(pop_mat[obs_time_inds, "I"] == 0)) {
                    # find the index of the first observation time where the rates are 0
                    end_ind <- which(obstimes > pop_mat[which(pop_mat[,"I"] == 0)[1], "time"])[1]
                    
                    # trim the observation times, observation matrix, and data matrix
                    obstimes <- obstimes[1:end_ind]
                    obs_mat <- obs_mat[1:end_ind, ]

                    # reset the final observation time in the configuration matrix
                    pop_mat <- pop_mat[pop_mat[,"time"] <= max(obstimes),]
          }
          
          return(list(pop_mat = pop_mat, obs_mat = obs_mat))
}