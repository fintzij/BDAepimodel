// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Update eigen values, vectors, and inverse matrices for irms
//'
//' @param tpm_prods array of transition probability matrix products
//' @param tpms array of transition probability matrices
//' @param pop_mat population level bookkeeping matrix
//' @param obstimes vector of observation times
//'
//' @return Updated eigenvalues, eigenvectors, and inverse matrices
// [[Rcpp::export]]
void tpmProdSeqs(Rcpp::NumericVector& tpm_prods, Rcpp::NumericVector& tpms, const Rcpp::IntegerVector obs_time_inds) {

        // Get indices
        int n_obs = obs_time_inds.size();
        Rcpp::IntegerVector tpmDims = tpms.attr("dim");

        // Ensure obs_time_inds starts at 0
        Rcpp::IntegerVector obs_inds = obs_time_inds - 1;

        // Set array pointers
        arma::cube prod_arr(tpm_prods.begin(), tpmDims[0], tpmDims[1], tpmDims[2], false);
        arma::cube tpm_arr(tpms.begin(), tpmDims[0], tpmDims[1], tpmDims[2], false);

        // Generate tpm products and store them in the appropriate spots in the tpm product array
        for(int j = 0; j < (n_obs - 1); ++j) {

                prod_arr.slice(obs_inds[j+1] - 1) = tpm_arr.slice(obs_inds[j+1] - 1);

                for(int k = (obs_inds[j+1] - 2); k > (obs_inds[j] - 1); --k) {

                        prod_arr.slice(k) = tpm_arr.slice(k) * prod_arr.slice(k + 1);
                }

        }
}