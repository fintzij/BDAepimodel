// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>

using namespace arma;
using namespace Rcpp;

//' Update eigen values, vectors, and inverse matrices for irms
//'
//' @param eigenvals matrix of eigenvalues with columns corresponding to irms
//' @param eigenvecs array of eigenvectors with third dimension corresponding to
//'   irms
//' @param inversevecs array of inverses of matrices of eigenvectors
//' @param irm_array array of rate matrices
//' @param complex_eigs logical vector indicating whether any eigen values are complex
//'
//' @return Updated eigenvalues, eigenvectors, and inverse matrices
//' @export
//' 
// [[Rcpp::export]]
void buildEigenArray(arma::mat& real_eigenvals, arma::mat& imag_eigenvals, arma::cube& eigenvecs, arma::cube& inversevecs, arma::cube& irm_array, Rcpp::IntegerVector& n_real_eigs) {

          // Get dimensions of the irm array
          int num_rows = irm_array.n_rows;
          int num_cols = irm_array.n_cols;
          int num_irms = irm_array.n_slices;
          int n_vals   = real_eigenvals.n_rows;
          
          // place-holder vector and matrix for complex eigenvalues/vectors
          arma::cx_vec valvec(n_vals);
          arma::cx_mat vecmat(num_rows, num_cols);

          // Compute new eigen decompositions
          for(int j = 0; j < num_irms; ++j) {

                    // Compute eigen decomposition and get inverse matrix
                    arma::eig_gen(valvec, vecmat, irm_array.slice(j));
                    
                    // count the number of real eigenvalues (modulo numeric error)
                    n_real_eigs[j] = arma::sum(abs(arma::imag(valvec)) < 0.000000000001);
          
                    // check whether any eigenvalues are complex
                    if(n_real_eigs[j] == n_vals) {
                              
                              // save the eigenvalues
                              real_eigenvals.col(j) = arma::real(valvec);
                              imag_eigenvals.col(j).zeros();

                              // save the eigenvectors and their inverse
                              eigenvecs.slice(j)   = arma::real(vecmat);
                              
                              // compute the pseudo-inverse of the eigenvector matrix
                              inversevecs.slice(j) = arma::pinv(eigenvecs.slice(j));
                              
                    } else {
                              // identify which eigenvalues are real and complex
                              // (0 is always an eigenvalue, so this is never NULL)
                              arma::uvec which_real    = arma::find(arma::imag(valvec) == 0.0);
                              arma::uvec which_complex = arma::find(arma::imag(valvec) != 0.0);
                              arma::uvec which_pos     = arma::find(arma::imag(valvec) > 0.0);
                              int n_complex            = which_pos.n_elem;

                              // save the real parts of the eigenvalues 
                              real_eigenvals(arma::span(0,n_real_eigs[j]-1), j)       = arma::real(valvec(which_real));
                              real_eigenvals(arma::span(n_real_eigs[j], n_vals-1), j) = arma::real(valvec(which_complex));
                              
                              // save the real eigenvectors
                              eigenvecs.slice(j).cols(0, n_real_eigs[j]-1) = arma::real(vecmat.cols(which_real));
                              
                              // save the complex eigenvalues
                              imag_eigenvals(arma::span(0,n_real_eigs[j]-1), j).zeros();
                              imag_eigenvals(arma::span(n_real_eigs[j], n_vals-1), j) = arma::imag(valvec(which_complex));
                              
                              // save the real and imaginary parts of positive conjugate eigenvector
                              for(int k = 0; k < n_complex; ++k) {
                                        eigenvecs.slice(j).col(n_real_eigs[j] + k)     = arma::real(vecmat.col(which_pos[k]));
                                        eigenvecs.slice(j).col(n_real_eigs[j] + k + 1) = arma::imag(vecmat.col(which_pos[k]));
                              }
                              
                              // compute the pseudo-inverse of the eigenvector matrix
                              inversevecs.slice(j) = arma::pinv(eigenvecs.slice(j));
                    }
          }
}
