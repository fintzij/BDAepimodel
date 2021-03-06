% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{simulateSIR}
\alias{simulateSIR}
\title{simulateSIR function, included in this file because of difficulty including sample.h in two separate files.
Rcpp code to simulate a general stochastic epidemic and binomial counts.}
\usage{
simulateSIR(obstimes, params, init_config)
}
\arguments{
\item{obstimes}{vector of observation times}

\item{params}{vector of parameters: beta, mu}

\item{init_config}{initial configuration of compartment counts}

\item{popsize}{population size}
}
\value{
updated bookeeping objects
}
\description{
simulateSIR function, included in this file because of difficulty including sample.h in two separate files.
Rcpp code to simulate a general stochastic epidemic and binomial counts.
}
