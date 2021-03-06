% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_tpm.R
\name{build_tpm}
\alias{build_tpm}
\title{Construct a transition probability matrix from}
\usage{
build_tpm(values, vectors, inv_vectors, t0, t1)
}
\arguments{
\item{values}{vector of eigenvalues}

\item{vectors}{matrix of eigenvectors}

\item{inv_vectors}{inverse of matrix of eigenvectors}

\item{t0}{left endpoint of the interval}

\item{t1}{right endpoint of the interval}
}
\value{
transition probability matrix
}
\description{
If (t1 - t0) is small, there could be numerical underflow leading to negative
entries in the tpm. If so, the negative values are set to 0.
}
