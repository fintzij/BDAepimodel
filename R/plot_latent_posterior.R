#' Plot the posterior distribution of the latent paths for measured states.
#' 
#' @param epimodel fitted epimodel object
#' @param states character vector of states for which to plot the posterior
#' @param times sequence of times at which to compute the posterior
#'   
#' @return plot of the posterior distribution of the latent paths
#' @export
#' 
plot_latent_posterior <- function(epimodel, states, times) {
          
          if(missing(states)) states <- epimodel$model$meas_vars
          if(missing(times)) times <- seq(min(epimodel$model$obstimes), max(epimodel$model$obstimes), length = 100)
          
          counts <- matrix(0, nrow = length(times) * length(states), ncol = length(epimodel$results$configs)) 
          
          for(s in 1:length(states)) {
                    for(k in 1:ncol(counts)) {
                              for(j in 1:length(times)) {
                                        
                                        config <- epimodel$results$configs[[k]]
                                        
                                        counts[j,k] <- config[sum(config[,"time"] <= times[j]), states[s]]
                              }
                    }
          }
          
          trajecs_dist <- data.frame(time = rep(times, length(states)), state = rep(states, length(times)), mean = 0, lower = 0, upper = 0)
          trajecs_dist$mean <- rowMeans(counts) 
          trajecs_se <- sqrt(sapply(apply(counts, 1, mcmc::initseq), "[[", "var.pos") / ncol(counts))
          trajecs_dist$lower <- trajecs_dist$mean - 1.96*trajecs_se
          trajecs_dist$upper <- trajecs_dist$mean + 1.96*trajecs_se
          
          return(print(ggplot2::ggplot(trajecs_dist, ggplot2::aes(x = time, y = mean, colour = state)) + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = state), alpha = 0.2) + ggplot2::geom_line() + ggplot2::theme_bw() + ggplot2::labs(title = "Posterior distribution of latent process", x = "Count", y = "Time")))
}