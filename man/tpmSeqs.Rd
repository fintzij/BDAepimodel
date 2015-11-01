% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{tpmSeqs}
\alias{tpmSeqs}
\title{Update eigen values, vectors, and inverse matrices for irms}
\usage{
tpmSeqs(tpms, pop_mat, eigen_vals, eigen_vecs, inverse_vecs, irm_keys)
}
\arguments{
\item{tpms}{array of transition probability matrices to be modified}

\item{pop_mat}{population level bookkeeping matrix}

\item{eigen_vals}{matrix of eigen values of irms}

\item{eigen_vecs}{array of eigen vectors of rate matrices}

\item{inverse_vecs}{array of inverse matrices of eigen vectors}

\item{irm_keys}{vector of irm array keys}
}
\value{
Updated eigenvalues, eigenvectors, and inverse matrices
}
\description{
Update eigen values, vectors, and inverse matrices for irms
}
