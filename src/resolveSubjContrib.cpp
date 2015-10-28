// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Add or remove the subject contribution from the compartment counts
//'
//' @param pop_mat population bookkeeping matrix in epimodel list
//' @param row_inds vector of row indices, 1:epimodel$ind_final_config
//' @param col_inds column inds given by the subject path in epimodel$config_mat
//' @param insertion Logical indicating whether the resolution should be the
//'     insertion or the removal of a path
//'
//' @return updated compartment counts
// [[Rcpp::export]]
void resolveSubjContrib(Rcpp::NumericVector& pop_mat, const Rcpp::IntegerVector row_inds, const Rcpp::IntegerVector col_inds, bool insertion) {

                Rcpp::IntegerVector pop_matDims = pop_mat.attr("dim");
                arma::mat newmat(pop_mat.begin(), pop_matDims[0], pop_matDims[1], false);
                int n = row_inds.size();

                if(insertion) {
                        for(int j = 0; j < n; ++j) {
                                newmat(row_inds[j] - 1, col_inds[j] +2) += 1; // there are always three columns before the states
                        }
                } else {
                        for(int j = 0; j < n; ++j) {
                                newmat(row_inds[j] - 1, col_inds[j] +2) -= 1; // there are always three columns before the states
                        }
                }

}