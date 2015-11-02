// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Get the irm keys for the compartment counts in a population level
//' bookkeeping matrix
//'
//' @param subject
//' @param pop_mat population level bookkeeping matrix
//' @param subj_path integer vector giving the subject path
//' @param irm_array array of rate matrices
//' @param keys vector of irm keys
//' @param loglik boolean indicating whether to return the likelihood or
//'     log-likelihood
//'
//' @return subject level likelihood or log-likelihood
// [[Rcpp::export]]
double subjectLikelihood(const int subject, const arma::mat& pop_mat, const arma::uvec& subj_path, Rcpp::NumericVector& irm_array, const Rcpp::NumericVector& initdist, const Rcpp::IntegerVector& keys, bool loglik) {

          // get number of intervals, and number of endpoints
          int final_ind = keys.size();
          arma::uword n_intervals = final_ind - 1;

          // Instatiate array pointers
          Rcpp::IntegerVector irmDims = irm_array.attr("dim");
          arma::cube irm(irm_array.begin(), irmDims[0], irmDims[1], irmDims[2], false);

          // initialize the likelihood
          double subj_lik = log(initdist[subj_path[0] - 1]);

          // compute the likelihood of the path
          for(uword j = 0; j < n_intervals; ++j) {
                    
                    if(subj_path[j] == subj_path[j+1]) {
                    
                              subj_lik += irm(subj_path[j] - 1, subj_path[j] - 1, keys[j] - 1) * (pop_mat(j + 1, 0) - pop_mat(j, 0));
                    
                    
                    } else {
                    
                              subj_lik += log(irm(subj_path[j] - 1, subj_path[j+1] - 1, keys[j] - 1)) + 
                                        irm(subj_path[j] - 1, subj_path[j] - 1, keys[j] - 1) * (pop_mat(j + 1, 0) - pop_mat(j, 0));
                    
                    }
          
          }

          if(loglik) {
                    return subj_lik;
          } else {
                    return exp(subj_lik);
          }
}