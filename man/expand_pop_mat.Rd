% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expand_pop_mat.R
\name{expand_pop_mat}
\alias{expand_pop_mat}
\title{Add buffer rows to a pop_mat matrix to enable in-place substitution of 
configurations.}
\usage{
expand_pop_mat(epimodel, buffer_size = NULL)
}
\arguments{
\item{buffer_size}{number of buffer rows to insert, defaulting to 30 percent 
of the number of rows in the current configuration matrix unless otherwise 
specified.}
}
\value{
new population level configuration matrix and index for the last row
  in the expanded configuration matrix (corresponding to the configuration at
  tmax)
}
\description{
Each buffer row is made up of NA entries that are inserted at the end of the 
configuration matrix. This function is designed to operate within an epimodel
environment and therefore makes multiple replacements within that 
environment.
}
