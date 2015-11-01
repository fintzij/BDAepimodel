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
arma::mat reorderMat(Rcpp::NumericMatrix& oldmtx, const arma::uvec& ord) {
          
          // Get matrix dimensions
          Rcpp::IntegerVector mtxDims = oldmtx.attr("dim");

          // shift the row indexing to start at the 0th element
          arma::uvec roword = ord - 1;
          
          // reorder the rows
          arma::mat newmtx(oldmtx.begin(), mtxDims[0], mtxDims[1], false);
          
          // return the matrix
          return newmtx.rows(roword);

}