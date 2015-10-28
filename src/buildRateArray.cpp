// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Update an array of rate matrices with the current rates
//'
//' @param irm_array array of rate matrices to be updated
//' @param rates matrix of rates
//' @param flow_inds matrix of indices in irm for where to place each rate
//'
//' @return updated array of rate matrices
// [[Rcpp::export]]
void buildRateArray(Rcpp::NumericVector& irm_array, const Rcpp::NumericMatrix& rates, const Rcpp::NumericMatrix& flow_inds) {

        // Get dimensions
        Rcpp::IntegerVector arrayDims = irm_array.attr("dim");
        int numRates = rates.ncol();

        // instatiate array pointers
        arma::cube arr(irm_array.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

        // update irm elements
        for(int j = 0; j < arrayDims[2]; ++j) {

                // set diagonal elements to 0
                arr.slice(j).diag() = arma::zeros(arrayDims[0]);

                // set off diagonal elements to the rates
                for(int k = 0; k < numRates; ++k) {
                        arr.slice(j)(flow_inds(k, 0) - 1, flow_inds(k, 1) - 1) = rates(j, k);
                }

                // set hazards to be the negative of the sum of off diagonals
                arr.slice(j) += arma::diagmat(-sum(arr.slice(j), 1));
        }
}

