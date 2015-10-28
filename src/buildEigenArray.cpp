// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Update eigen values, vectors, and inverse matrices for irms
//'
//' @param eigenvals matrix of eigenvalues with columns corresponding to irms
//' @param eigenvecs array of eigenvectors with third dimension corresponding to
//'   irms
//' @param inversevecs array of inverses of matrices of eigenvectors
//' @param irm_array array of rate matrices
//'
//' @return Updated eigenvalues, eigenvectors, and inverse matrices
// [[Rcpp::export]]
void buildEigenArray(arma::mat& eigenvals, Rcpp::NumericVector& eigenvecs, Rcpp::NumericVector& inversevecs, Rcpp::NumericVector& irm_array) {

        // Get dimensions of the irm array
        Rcpp::IntegerVector irmDims = irm_array.attr("dim"); // the vector cubes also have the same dimensions
        int n_vals = eigenvals.n_rows;

        // instatiate pointers
        arma::cube vec_arr(eigenvecs.begin(), irmDims[0], irmDims[1], irmDims[2], false);
        arma::cube inv_arr(inversevecs.begin(), irmDims[0], irmDims[1], irmDims[2], false);
        arma::cube irm_arr(irm_array.begin(), irmDims[0], irmDims[1], irmDims[2], false);

        // place-holder vector and matrix for complex eigenvalues/vectors
        arma::cx_vec valvec(n_vals);
        arma::cx_mat vecmat(irmDims[0], irmDims[1]);
        arma::mat invmat(irmDims[0], irmDims[1]);

        // Compute new eigen decompositions
        for(int j = 0; j < irmDims[2]; ++j) {

                // Compute eigen decomposition and get inverse matrix
                arma::eig_gen(valvec, vecmat, irm_arr.slice(j));
                invmat = arma::pinv(arma::real(vecmat));

                // Make replacements
                eigenvals.col(j) = arma::real(valvec);
                vec_arr.slice(j) = arma::real(vecmat);
                inv_arr.slice(j) = invmat;
        }
}