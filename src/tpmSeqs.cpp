// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Update eigen values, vectors, and inverse matrices for irms
//'
//' @param tpms array of transition probability matrices to be modified
//' @param pop_mat population level bookkeeping matrix
//' @param eigen_vals matrix of eigen values of irms
//' @param eigen_vecs array of eigen vectors of rate matrices
//' @param inverse_vecs array of inverse matrices of eigen vectors
//' @param irm_keys vector of irm array keys
//'
//' @return Updated eigenvalues, eigenvectors, and inverse matrices
// [[Rcpp::export]]
void tpmSeqs(Rcpp::NumericVector& tpms, arma::mat& pop_mat, arma::mat eigen_vals, Rcpp::NumericVector& eigen_vecs, Rcpp::NumericVector& inverse_vecs, Rcpp::IntegerVector& irm_keys) {

        // Set constants
        int n_intervals = irm_keys.size() - 1;
        Rcpp::IntegerVector arrayDims = eigen_vecs.attr("dim");
        Rcpp::IntegerVector tpmDims = tpms.attr("dim");

        // Set keys to start at 0
        Rcpp::IntegerVector keys = irm_keys - 1;

        // Set pointers
        arma::cube tpm_arr(tpms.begin(), tpmDims[0], tpmDims[1], tpmDims[2], false);
        arma::cube vec_arr(eigen_vecs.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
        arma::cube inv_arr(inverse_vecs.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

        // get time diffs
        arma::vec times_all = pop_mat.col(0);
        arma::vec timediffs = diff(times_all.subvec(0, n_intervals));

        // compute new tpms
        for(int j = 0; j < n_intervals; ++j) {

                tpm_arr.slice(j) = vec_arr.slice(keys[j]) * (arma::diagmat(exp(timediffs[j] * eigen_vals.col(keys[j]))) * inv_arr.slice(keys[j]));

                tpm_arr.slice(j).elem(find(tpm_arr.slice(j) < 1e-13)).zeros(); // set negatives and numerical errors to zero

        }
}