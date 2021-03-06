% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_params.R
\name{update_params}
\alias{update_params}
\title{Replace the vector of parameters in an epimodel object with new values and 
update the population level likelihood, the measurement process 
likelihood, and the eigen decomposition objects for irms.}
\usage{
update_params(epimodel, params, pop_likelihood = NULL,
  obs_likelihood = NULL)
}
\arguments{
\item{obs_likelihood}{}
}
\value{
updated parameter vector and likelihood objects in the epimodel
  environment
}
\description{
Replace the vector of parameters in an epimodel object with new values and 
update the population level likelihood, the measurement process 
likelihood, and the eigen decomposition objects for irms.
}
