// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Construct a transition probability matrix for a rate matrix with complex 
//' eigenvalues using the real canonical form of the decomposition.
//' 
//' @param vals vector of eigenvalues of type cx_vec
//' @param vecs matrix of eigenvectors (positive complex conjugate)
//' @param inv_vecs inverse matrix of eigenvectors
//' @param dt time interval
//'
//' @return Updated eigenvalues, eigenvectors, and inverse matrices
// [[Rcpp::export]]
arma::mat complexTPM(const arma::vec& real_vals, const arma::vec& imag_vals, const arma::mat& vecs, const arma::mat& inv_vecs, double dt, int n_real) {
          
          // get the dimensions of the tpm (n_vals x n_vals)
          int n_vals  = real_vals.n_elem;
          int n_pairs = (n_vals - n_real)/2;

          // initialize the exponential matrix and a submatrix for the complex matrix exponentials
          arma::mat exp_Bt(n_vals, n_vals, arma::fill::zeros);
          arma::mat exp_sub(2,2, arma::fill::zeros);
          double at = 0.0;
          double bt = 0.0;

          // fill out the diagonal terms of exp_Bt for the real eigenvalues
          for(int k=0; k < n_real; ++k) {
                    exp_Bt(k,k) = exp(real_vals[k] * dt);
          }
          
          // fill out the block diagonal terms of exp_Bt for the complex eigenvalues
          for(int k=0; k < n_pairs; ++k) {
                    at = real_vals[2*k + n_real] * dt;
                    bt = imag_vals[2*k + n_real] * dt;
                    exp_sub(0,0) = cos(bt);
                    exp_sub(0,1) = sin(bt);
                    exp_sub(1,0) = -sin(bt);
                    exp_sub(1,1) = cos(bt);
                    exp_Bt(arma::span(n_real+k, n_real+k+1), arma::span(n_real+k, n_real+k+1)) = exp(at) * exp_sub;
          }
          
          // return e^{Qt} = P * exp(Bt) * P^-1
          return (vecs * (exp_Bt * inv_vecs));
}

//' Update eigen values, vectors, and inverse matrices for irms
//'
//' @param tpms array of transition probability matrices to be modified
//' @param pop_mat population level bookkeeping matrix
//' @param real_eigen_vals matrix with the real parts of the eigen values
//' @param imag_eigen_vals matrix with the imaginary parts of the eigenvalues
//' @param eigen_vecs array of eigen vectors of rate matrices
//' @param inverse_vecs array of inverse matrices of eigen vectors
//' @param irm_keys vector of irm array keys
//' @param n_real_eigs integer vector indicating how many eigen values are real
//'
//' @return Updated eigenvalues, eigenvectors, and inverse matrices
// [[Rcpp::export]]
void tpmSeqs(arma::cube& tpms, const arma::mat& pop_mat, const arma::mat& real_eigen_vals, const arma::mat& imag_eigen_vals, const arma::cube& eigen_vecs, const arma::cube& inverse_vecs, const Rcpp::IntegerVector& irm_keys, const Rcpp::IntegerVector& n_real_eigs, const arma::cube& irms) {

          // Set constants
          int n_intervals = irm_keys.size() - 1;
          int n_vals      = real_eigen_vals.n_rows;
          
          // Set keys to start at 0
          Rcpp::IntegerVector keys = irm_keys - 1;
          
          // get time diffs
          arma::vec times_all = pop_mat.col(0);
          arma::vec timediffs = diff(times_all.subvec(0, n_intervals));
          
          // compute new tpms
          for(int j = 0; j < n_intervals; ++j) {
                  
                    if(n_real_eigs[keys[j]] == n_vals) {
                              tpms.slice(j) = eigen_vecs.slice(keys[j]) * (arma::diagmat(exp(timediffs[j] * arma::real(real_eigen_vals.col(keys[j])))) * inverse_vecs.slice(keys[j]));
                    } else {
                              tpms.slice(j) = complexTPM(real_eigen_vals.col(keys[j]), imag_eigen_vals.col(keys[j]), eigen_vecs.slice(keys[j]), inverse_vecs.slice(keys[j]), timediffs[j], n_real_eigs[keys[j]]);
                    }
                    
                    // set negatives and numerical errors to zero
                    tpms.slice(j).elem(find(tpms.slice(j) < 1e-13)).zeros(); 
          }
}