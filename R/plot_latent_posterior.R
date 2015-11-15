#' Plot the posterior distribution of the latent paths for measured states.
#' 
#' @param epimodel fitted epimodel object
#' @param states character vector of states for which to plot the posterior
#' @param times sequence of times at which to compute the posterior
#' @param proportion if true, plot the proportion instead of the counts
#' @param truth true latent trajectory
#' @param true_popsize true population size, defaults to the population size
#'   specified in the epimodel object unless otherwise provided.
#'   
#' @return plot of the posterior distribution of the latent paths
#' @export
#' 
plot_latent_posterior <- function(epimodel, states, times, proportion = FALSE, truth = NULL, true_popsize = NULL) {
          
          if(missing(states)) states <- epimodel$model$meas_vars
          if(missing(times)) times <- seq(min(epimodel$model$obstimes), max(epimodel$model$obstimes), length = 201)
          
          counts <- matrix(0, nrow = length(times) * length(states), ncol = length(epimodel$results$configs)) 
          
          for(s in 1:length(states)) {
                    for(k in 1:ncol(counts)) {
                              for(j in 1:length(times)) {
                                        
                                        config <- epimodel$results$configs[[k]]
                                        
                                        counts[j + length(times) * (s - 1), k] <- config[sum(config[,"time"] <= times[j]), states[s]]
                              }
                    }
          }
          
          if(proportion) {
                    counts <- counts / epimodel$model$popsize
          }
          
          trajecs_dist <- data.frame(time = rep(times, length(states)), state = rep(states, length(times)), mean = 0, lower = 0, upper = 0)
          trajecs_dist$mean <- rowMeans(counts) 
          trajecs_se <- sqrt(sapply(apply(counts, 1, mcmc::initseq), "[[", "var.pos"))
          trajecs_dist$lower <- trajecs_dist$mean - 1.96*trajecs_se
          trajecs_dist$upper <- trajecs_dist$mean + 1.96*trajecs_se
          
          if(proportion) {
                    trajecs_post <- ggplot2::ggplot(trajecs_dist, ggplot2::aes(x = trajecs_dist$time, y = trajecs_dist$mean, colour = trajecs_dist$state)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = trajecs_dist$lower, ymax = trajecs_dist$upper, fill = trajecs_dist$state), alpha = 0.2) + ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::labs(title = "Posterior distribution of latent process", x = "Time", y = "Proportion") 
          } else {
                    trajecs_post <- ggplot2::ggplot(trajecs_dist, ggplot2::aes(x = trajecs_dist$time, y = trajecs_dist$mean, colour = trajecs_dist$state)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = trajecs_dist$lower, ymax = trajecs_dist$upper, fill = trajecs_dist$state), alpha = 0.2) + ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::labs(title = "Posterior distribution of latent process", x = "Time", y = "Count") 
                    
          }

          if(!is.null(truth)) {
                    if(is.null(true_popsize)) true_popsize <- epimodel$model$popsize
                    
                    true_count = data.frame(time = rep(times, length(states)), state = rep(states, length(times)), count = 0)
                    
                    for(s in 1:length(states)) {
                              
                              for(j in 1:length(times)) {
                                        
                                       true_count[j + length(times) * (s - 1), "count"]  <- truth[sum(truth[,"time"] <= times[j]), states[s]]
                              }
                    }
                    
                    if(proportion) true_count$count  <- true_count$count / true_popsize
                    
                    trajecs_post <- trajecs_post + ggplot2::geom_line(data = true_count, aes(x = time, y = count), colour = "black") + guides(gill = guide_legend(title = "State"))
                    
          } else {
                    trajecs_post <- trajecs_post + guides(fill = guide_legend(title ="State"))
          }
          
          return(print(ggplot2::ggplot(trajecs_dist, ggplot2::aes(x = trajecs_dist$time, y = trajecs_dist$mean, colour = trajecs_dist$state)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = trajecs_dist$lower, ymax = trajecs_dist$upper, fill = trajecs_dist$state), alpha = 0.2) + ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::labs(title = "Posterior distribution of latent process", x = "Count", y = "Time")))
}