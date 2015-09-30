% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/check_irm.R
\name{check_irm}
\alias{check_irm}
\title{Checks that all of the IRMs that are required exist in the .irm environment
by checking keys against the .irm object list.}
\usage{
check_irm(epimodel)
}
\arguments{
\item{epimodel}{list of bookkeeping and model objects, or an analagous
environment. If not supplied, the list will be generated. The epimodel list
is coerced into an environment, .epimodel, for internal use.}
}
\value{
Instatiated rate matrix in the .irm environment and its eigen
  decomposition in the .eigen environment
}
\description{
Checks that all of the IRMs that are required exist in the .irm environment
by checking keys against the .irm object list.
}
