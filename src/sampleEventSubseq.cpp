// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

//' Update eigen values, vectors, and inverse matrices for irms
//'
//' @param tpms array of transition probability matrices to be modified
//' @param pop_mat population level bookkeeping matrix
//' @param eigen_vals matrix of eigen values of irms
//' @param eigen_vecs array of eigen vectors of rate matrices
//' @param inverse_vecs array of inverse matrices of eigen vectors
//' @param irm_keys vector of irm array keys
//'
//' @return Updated eigenvalues, eigenvectors, and inverse matrices
// [[Rcpp::export]]
Rcpp::IntegerVector sampleEventSubseq(Rcpp::IntegerVector& path, Rcpp::NumericVector& tpms, Rcpp::NumericVector& tpm_prods, const int init_ind, const int final_ind) {
          
          // set the initial and final inds
          int ind_start = init_ind - 1;
          int ind_end = final_ind -1;
          
          // set the initial and final states
          int state_cur = path[ind_start];
          int state_end = path[ind_end];
          Rcpp::IntegerVector next_state(1);
          
          // get tpm dimensions
          Rcpp::IntegerVector tpmDims = tpms.attr("dim");
          
          // set tpm pointers and path pointer
          arma::cube tpm_arr(tpms.begin(), tpmDims[0], tpmDims[1], tpmDims[2], false);
          arma::cube prod_arr(tpm_prods.begin(), tpmDims[0], tpmDims[1], tpmDims[2], false);
          
          // create vector for state probabilities
          Rcpp::NumericVector state_probs(tpmDims[0]);
          Rcpp::IntegerVector states = Rcpp::seq_len(tpmDims[0]);
          
          // Sample the status at observation times
          for(int j = (ind_start+1); j < ind_end; ++j) {
                    
                    state_probs = (tpm_arr.slice(j - 1).row(state_cur - 1).t() % prod_arr.slice(j).col(state_end - 1)) / prod_arr.slice(j - 1).at(state_cur - 1, state_end - 1);
                    
                    next_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs);
                    
                    path[j] = state_cur = int(next_state[0]);
          }
          
          return path;
}