// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Add or remove the subject contribution from the compartment counts
//'
//' @param pop_mat population bookkeeping matrix in epimodel list
//' @param config_mat subject level configuration bookkeeping matrix
//' @param subj_inds indices for rows relating to the subject's events
//'
//' @return updated compartment counts
// [[Rcpp::export]]
void removeSubjPath(Rcpp::NumericVector& pop_mat, Rcpp::NumericVector& config_mat, Rcpp::IntegerVector subj_inds) {


}