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
          
          # if which_compartments is not NULL, ensure trajec is TRUE
          if(!is.null(which_compartments)) {
                    config <- TRUE
          } else if(is.null(which_compartments) & config) {
                    which_compartments <- epimodel$states
          } else if(is.null(which_compartments) & !config) {
                    which_compartments <- NULL
          }
          
          if(obs & config) {
                    config_melt <- reshape2::melt(data.frame(epimodel$pop_mat[complete.cases(epimodel$pop_mat),c("time", which_compartments)]), id.vars = "time")
                    obs_melt <- reshape2::melt(as.data.frame(epimodel$dat[,c("time", epimodel$meas_vars)]), id.vars = epimodel$time_var)
                    
                    epiplot <- ggplot2::ggplot(config_melt, ggplot2::aes(x = config_melt$time, y = config_melt$value, colour = config_melt$variable)) + ggplot2::geom_step() + ggplot2::theme_bw() + ggplot2::labs(x = "Time", y = "Count") + ggplot2::scale_colour_discrete(guide = ggplot2::guide_legend(title = "True compartment count")) + ggplot2::geom_point(data = obs_melt, ggplot2::aes(x = obs_melt$time, y = obs_melt$value, colour = obs_melt$variable, shape = obs_melt$variable)) + ggplot2::scale_shape_discrete(guide = ggplot2::guide_legend(title = "Observed count"))
          }
          
          if(obs & !config) {
                    obs_melt <- reshape2::melt(as.data.frame(epimodel$dat[,c("time", epimodel$meas_vars)]), id.vars = epimodel$time_var)
                    
                    epiplot <- ggplot2::ggplot(obs_melt, ggplot2::aes(x = obs_melt$time, y = obs_melt$value, colour = obs_melt$variable)) + ggplot2::geom_point() + ggplot2::theme_bw() + ggplot2::labs(x = "Time", y = "Count")
          }
          
          if(!obs & config) {
                    config_melt <- reshape2::melt(as.data.frame(epimodel$pop_mat[,c("time", which_compartments)]), id.vars = "time")
                    
                    epiplot <- ggplot2::ggplot(config_melt, ggplot2::aes(x = config_melt$time, y = config_melt$value, colour = config_melt$variable)) + ggplot2::geom_step() + ggplot2::theme_bw() + ggplot2::labs(x = "Time", y = "Count") + ggplot2::scale_colour_discrete(guide = ggplot2::guide_legend(title = "Compartment"))
          }
          
          return(print(epiplot))
}