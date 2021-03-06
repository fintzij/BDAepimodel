% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{complexTPM}
\alias{complexTPM}
\title{Construct a transition probability matrix for a rate matrix with complex 
eigenvalues using the real canonical form of the decomposition.}
\usage{
complexTPM(real_vals, imag_vals, vecs, inv_vecs, dt, n_real)
}
\arguments{
\item{vecs}{matrix of eigenvectors (positive complex conjugate)}

\item{inv_vecs}{inverse matrix of eigenvectors}

\item{dt}{time interval}

\item{vals}{vector of eigenvalues of type cx_vec}
}
\value{
Updated eigenvalues, eigenvectors, and inverse matrices
}
\description{
Construct a transition probability matrix for a rate matrix with complex 
eigenvalues using the real canonical form of the decomposition.
}
