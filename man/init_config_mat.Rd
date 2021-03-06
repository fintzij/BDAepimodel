% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_config_mat.R
\name{init_config_mat}
\alias{init_config_mat}
\title{initializes the bookkeeping object for population level counts}
\usage{
init_config_mat(epimodel, init_state, t0, tmax)
}
\arguments{
\item{epimodel}{an epimodel list}

\item{t0}{the first observation time}

\item{tmax}{the final observation time}
}
\value{
initialized matrix with columns to store event times and counts of
  individuals in each compartment. The first and last rows of the matrix
  always correspond to time 0 and tmax. The state at tmax is initialized to
  init_state.
}
\description{
initializes the bookkeeping object for population level counts
}
\examples{
init_state <- c(S = 45, I = 5, R = 0)
tmax <- 5
init_pop_traj(init_state, tmax)

}
