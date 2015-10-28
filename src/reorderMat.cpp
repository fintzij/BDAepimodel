// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Reorder the rows of a matrix
//'
//' @param mtx matrix to be reordered
//' @param ord row order
//'
//' @return new array joining cubes
// [[Rcpp::export]]
arma::mat reorderMat(arma::mat& oldmtx, arma::uvec& ord) {

        // shift the row indexing to start at the 0th element
        arma::uvec roword = ord - 1;

        // reorder the rows
        arma::mat newmtx = oldmtx.rows(roword);

        // return the matrix
        return newmtx;

}