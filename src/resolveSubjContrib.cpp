// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Add or remove the subject contribution from the compartment counts
//'
//' @param pop_mat population bookkeeping matrix in epimodel list
//' @param row_inds vector of row indices, 1:epimodel$ind_final_config
//' @param subj_path column inds given by the subject path
//' @param insertion Logical indicating whether the resolution should be the
//'     insertion or the removal of a path
//'
//' @return updated compartment counts
// [[Rcpp::export]]
void resolveSubjContrib(Rcpp::NumericMatrix& pop_mat, const int ind_final_config, const Rcpp::IntegerVector& subj_path, bool insertion) {

                if(insertion) {
                        for(int j = 0; j < ind_final_config; ++j) {
                                pop_mat(j, subj_path[j] +2) += 1; // there are always three columns before the states
                        }
                } else {
                        for(int j = 0; j < ind_final_config; ++j) {
                                pop_mat(j, subj_path[j] +2) -= 1; // there are always three columns before the states
                        }
                }
}