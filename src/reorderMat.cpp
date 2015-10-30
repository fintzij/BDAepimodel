// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Reorder the rows of a matrix
//'
//' @param mtx matrix to be reordered
//' @param ord row order
//'
//' @return matrix with rows permuted according to ord
// [[Rcpp::export]]
arma::mat reorderMat(const arma::mat& oldmtx, const arma::uvec& ord) {

        // shift the row indexing to start at the 0th element
        arma::uvec roword = ord - 1;

        // reorder the rows
        arma::mat newmtx = oldmtx.rows(roword);

        // return the matrix
        return newmtx;

}