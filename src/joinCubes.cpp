// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Join two armadillo cubes
//'
//' @param cube1
//' @param cube2
//'
//' @return new array joining cubes
// [[Rcpp::export]]
arma::cube joinCubes(Rcpp::NumericVector& firstcube, Rcpp::NumericVector& secondcube) {

        // Get dimensions
        Rcpp::IntegerVector firstDims = firstcube.attr("dim");
        Rcpp::IntegerVector secondDims = secondcube.attr("dim");

        // Instatiate pointers
        arma::cube firstarr(firstcube.begin(), firstDims[0], firstDims[1], firstDims[2], false);
        arma::cube secondarr(secondcube.begin(), secondDims[0], secondDims[1], secondDims[2], false);

        // join the cubes
        arma::cube newcube = arma::join_slices(firstarr, secondarr);

        return newcube;
}