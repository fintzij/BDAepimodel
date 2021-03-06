% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{buildRateArray}
\alias{buildRateArray}
\title{Update an array of rate matrices with the current rates}
\usage{
buildRateArray(irm_array, rates, flow_inds)
}
\arguments{
\item{irm_array}{array of rate matrices to be updated}

\item{rates}{matrix of rates}

\item{flow_inds}{matrix of indices in irm for where to place each rate}
}
\value{
updated array of rate matrices
}
\description{
Update an array of rate matrices with the current rates
}
