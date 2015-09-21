#' Plot the current trajectory and data in an epimodel object (uses ggplot2)
#'
#' @inheritParams simulate_epimodel
#' @param obs TRUE/FALSE for whether to plot the observed data 
#' @param config TRUE/FALSE for whether to plot the population level trajectory
#' @param which_compartments character vector specifying which compartments to plot (defaults to all compartments)
#'
#' @return plot of trajectories and/or data
#' @export
#'
plot_epimodel <- function(epimodel, obs = TRUE, config = TRUE, which_compartments = NULL) {
          
          require(ggplot2)
          require(reshape2)
          
          # if which_compartments is not NULL, ensure trajec is TRUE
          if(!is.null(which_compartments)) {
                    config <- TRUE
          } else if(is.null(which_compartments) & config) {
                    which_compartments <- epimodel$states
          } else if(is.null(which_compartments) & !config) {
                    which_compartments <- NULL
          }
          
          if(obs & config) {
                    config_melt <- melt(as.data.frame(epimodel$config_mat[complete.cases(epimodel$config_mat),c("time", which_compartments)]), id.vars = "time")
                    obs_melt <- melt(as.data.frame(epimodel$dat[,c("time", epimodel$meas_vars)]), id.vars = epimodel$time_var)
                    
                    epiplot <- ggplot(config_melt, aes(x = time, y = value, colour = variable)) + geom_step() + theme_bw() + labs(x = "Time", y = "Count") + scale_colour_discrete(guide = guide_legend(title = "True compartment count")) + geom_point(data = obs_melt, aes(x = time, y = value, colour = variable, shape = variable)) + scale_shape_discrete(guide = guide_legend(title = "Observed count"))
          }
          
          if(obs & !config) {
                    obs_melt <- melt(as.data.frame(epimodel$dat[,c("time", epimodel$meas_vars)]), id.vars = epimodel$time_var)
                    
                    epiplot <- ggplot(obs_melt, aes(x = time, y = value, colour = variable)) + geom_point() + theme_bw() + labs(x = "Time", y = "Count")
          }
          
          if(!obs & config) {
                    config_melt <- melt(as.data.frame(epimodel$config_mat[,c("time", which_compartments)]), id.vars = "time")
                    
                    epiplot <- ggplot(config_melt, aes(x = time, y = value, colour = variable)) + geom_step() + theme_bw() + labs(x = "Time", y = "Count") + scale_colour_discrete(guide = guide_legend(title = "Compartment"))
          }
          
          return(print(epiplot))
}