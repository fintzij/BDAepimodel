// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Construct forward-backward matrices
//'
//' @param fb_mats array of FB matrices to be updated
//' @param tpm_prods array of tpm products
//' @param emit_mat matrix of emission probabilities
//' @param initdist vector of initial state probabilities
//' @param obs_time_inds vector of observation time indices
//'
//' @return Updated array of FB matrices
// [[Rcpp::export]]

void buildFBMats(Rcpp::NumericVector& fb_mats, Rcpp::NumericVector& tpm_prods, arma::mat& emit_mat, arma::vec& initdist, Rcpp::IntegerVector& obs_time_inds) {
          
          // get dimensions of objects
          Rcpp::IntegerVector fbDims = fb_mats.attr("dim");
          Rcpp::IntegerVector prodDims = tpm_prods.attr("dim");

          // create pointers
          arma::cube fb_arr(fb_mats.begin(), fbDims[0], fbDims[1], fbDims[2], false);
          arma::cube prod_arr(tpm_prods.begin(), prodDims[0], prodDims[1], prodDims[2], false);
          
          // set observation time indices to start at 0
          Rcpp::IntegerVector obs_inds = obs_time_inds - 1;
          
          // Create the first FB matrix
          arma::vec pi_0 = arma::normalise(initdist % emit_mat.col(0), 1);
          fb_arr.slice(0) = (pi_0 * emit_mat.col(1).t()) % prod_arr.slice(0);
          fb_arr.slice(0) /= accu(fb_arr.slice(0));
          
          for(int j = 1; j < fbDims[2]; ++j) {
                    fb_arr.slice(j) = (sum(fb_arr.slice(j-1)).t() * emit_mat.col(j+1).t()) % prod_arr.slice(obs_inds[j]);
                    fb_arr.slice(j) /= accu(fb_arr.slice(j));
          }
}