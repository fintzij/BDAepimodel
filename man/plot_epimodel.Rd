% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_epimodel.R
\name{plot_epimodel}
\alias{plot_epimodel}
\title{Plot the current trajectory and data in an epimodel object (uses ggplot2)}
\usage{
plot_epimodel(epimodel, obs = TRUE, config = TRUE,
  which_compartments = NULL)
}
\arguments{
\item{obs}{TRUE/FALSE for whether to plot the observed data}

\item{config}{TRUE/FALSE for whether to plot the population level trajectory}

\item{which_compartments}{character vector specifying which compartments to plot (defaults to all compartments)}
}
\value{
plot of trajectories and/or data
}
\description{
Plot the current trajectory and data in an epimodel object (uses ggplot2)
}
