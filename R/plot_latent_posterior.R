#' Plot the posterior distribution of the latent paths for measured states.
#' 
#' @param epimodel fitted epimodel object
#' @param states character vector of states for which to plot the posterior
#' @param times sequence of times at which to compute the posterior
#' @param cm "med" if posterior median, 2.5 and 97.5 quantiles are desired,
#'   otherwise "mn" if posterior mean, and 95 percent CI are desired
#' @param proportion if true, plot the proportion instead of the counts
#' @param truth true latent trajectory
#' @param true_popsize true population size, defaults to the population size 
#'   specified in the epimodel object unless otherwise provided.
#' @param dat data matrix in the epimodel object
#'   
#' @return plot of the posterior distribution of the latent paths
#' @export 
plot_latent_posterior <- function(epimodel, states, times, cm = "med", proportion = FALSE, truth = NULL, true_popsize = NULL, dat = NULL) {
          
          if(missing(states)) states <- epimodel$meas_vars
          if(missing(times)) times <- seq(min(epimodel$obstimes), max(epimodel$obstimes), length = 201)
          
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
                    counts <- counts / epimodel$popsize
          }
          
          trajecs_dist <- data.frame(time = rep(times, length(states)), state = rep(states, length(times)), cm = 0, lower = 0, upper = 0)
          
          if(cm == "med") {
                    trajecs_dist$cm <- apply(counts, 1, median, na.rm = T)
                    
          } else if(cm == "mn") {
                    trajecs_dist$cm <- rowMeans(counts)
          }
          
          trajecs_dist$lower <- apply(counts, 1, quantile, 0.025)
          trajecs_dist$upper <- apply(counts, 1, quantile, 0.975)

          if(!is.null(truth)) {
                    true_count = data.frame(time = rep(times, length(states)), state = rep(states, length(times)), count = 0)
                    
                    for(s in 1:length(states)) {
                              
                              for(j in 1:length(times)) {
                                        
                                       true_count[j + length(times) * (s - 1), "count"]  <- truth[sum(truth[,"time"] <= times[j]), states[s]]
                              }
                    }
                    
                    if(proportion) true_count$count  <- true_count$count / true_popsize
                    
                    
          } 
          
          if(proportion & is.null(truth)) {
                    trajecs_post <- ggplot2::ggplot(trajecs_dist, ggplot2::aes(x = time, y = cm, colour = paste(state, "Augmented"))) + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = state), alpha = 0.2) + ggplot2::geom_line(size = 1)  + ggplot2::theme_bw() + ggplot2::scale_fill_discrete(guide="none") + ggplot2::labs(title = "Posterior distribution of latent process", x = "Time", y = "Proportion") + scale_color_discrete(guide = ggplot2::guide_legend(title = "State", override.aes = list(fill ="white"))) + ggplot2::theme(text = ggplot2::element_text(size = 15))
                    
          } else if(!proportion & is.null(truth)){
                    
                    trajecs_post <- ggplot2::ggplot(trajecs_dist, ggplot2::aes(x = time, y = cm, colour = paste(state, "Augmented"))) + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = state), alpha = 0.2) + ggplot2::geom_line(size = 1) + ggplot2::theme_bw() + ggplot2::scale_fill_discrete(guide="none") + ggplot2::labs(title = "Posterior distribution of latent process", x = "Time", y = "Count") + scale_color_discrete(guide = ggplot2::guide_legend(title = "State", override.aes = list(fill ="white"))) + ggplot2::theme(text = ggplot2::element_text(size = 15))
                    
          } else if(proportion & !is.null(truth)) {
                    trajecs_post <- ggplot2::ggplot(trajecs_dist, ggplot2::aes(x = time, y = cm, colour = paste(state, "Augmented"))) + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = state), alpha = 0.2) + ggplot2::geom_line(size = 1) + ggplot2::geom_step(data = true_count, ggplot2::aes(x = time, y = count, colour = paste(state,"Truth")), size = 1) + ggplot2::theme_bw() + ggplot2::scale_fill_discrete(guide="none") + ggplot2::labs(title = "Posterior distribution of latent process", x = "Time", y = "Proportion") + scale_color_discrete(guide = ggplot2::guide_legend(title = "State", override.aes = list(fill ="white"))) + ggplot2::theme(text = ggplot2::element_text(size = 15))
                    
          } else {
                    trajecs_post <- ggplot2::ggplot(trajecs_dist, ggplot2::aes(x = time, y = cm, colour = paste(state, "Augmented"))) + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = state), alpha = 0.2) + ggplot2::geom_line(size = 1) + ggplot2::geom_step(data = true_count, ggplot2::aes(x = time, y = count, colour = paste(state,"Truth")), size = 1) + ggplot2::theme_bw() + ggplot2::scale_fill_discrete(guide="none") + ggplot2::labs(title = "Posterior distribution of latent process", x = "Time", y = "Count") + scale_color_discrete(guide = ggplot2::guide_legend(title = "State", override.aes = list(fill ="white"))) + ggplot2::theme(text = ggplot2::element_text(size = 15))
          }
          
          if(!is.null(dat)) {
                    if(!proportion) {
                              trajecs_post <- trajecs_post + ggplot2::geom_point(aes(x = time, y = eval(parse(text = epimodel$meas_vars)), colour = paste(epimodel$model$meas_vars,"Observed Count")), size = 5, data = data.frame(epimodel$dat)) 
                    } else {
                              prop <- dat
                              prop[,epimodel$meas_vars] <- prop[,epimodel$meas_vars] / epimodel$popsize * 2
                              trajecs_post <- trajecs_post + ggplot2::geom_point(aes(x = time, y = eval(parse(text = epimodel$meas_vars)), colour = paste(epimodel$meas_vars,"Observed Count")), size = 5, data = data.frame(prop)) 
                    }
                    
          }
          
          return(print(trajecs_post))
}