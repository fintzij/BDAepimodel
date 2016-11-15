// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>

using namespace arma;
using namespace Rcpp;

//' Update eigen values, vectors, and inverse matrices analytically for the SIR 
//' model
//'
//' @param eigenvals matrix of eigenvalues with columns corresponding to irms
//' @param eigenvecs array of eigenvectors with third dimension corresponding to
//'   irms
//' @param inversevecs array of inverses of matrices of eigenvectors
//' @param irm_array array of rate matrices
//' @param complex_eigs logical vector indicating whether any eigen values are complex
//' @param first_calc boolean indicating whether this is the first time that the 
//'   eigen decomposition is being computed
//'
//' @return Updated eigenvalues, eigenvectors, and inverse matrices
//' @export
//' 
// [[Rcpp::export]]
void buildEigenArray_SIR(arma::mat& real_eigenvals, arma::mat& imag_eigenvals, arma::cube& eigenvecs, arma::cube& inversevecs, arma::cube& irm_array, Rcpp::IntegerVector& n_real_eigs, bool initial_calc) {
          
          // Get dimensions of the irm array
          int num_irms = irm_array.n_slices;

          // third eigenvalue is always 0 and third eigenvector is always 1
          if(initial_calc) {
                    real_eigenvals.zeros();
                    imag_eigenvals.zeros();
                    eigenvecs.each_slice() = arma::eye(3,3);    // set each each eigenvector matrix to diag(1,3)
                    eigenvecs.tube(0,2,arma::size(2,1)).ones(); // set the last eigenvector to 1
                    inversevecs.each_slice() = arma::eye(3,3);  // set each inverse matrix to diag(1,3)
                    inversevecs.tube(1,2).fill(-1);             // set the 2,3 element of the inverse vectors to -1
                    n_real_eigs.fill(3);                        // all three eigenvalues are real
          }
          
          // initialize parameter values
          double beta = 0;
          double mu = 0;
          
          // Compute new eigen decompositions
          for(int j = 0; j < num_irms; ++j) {
                    
                    beta = irm_array.slice(j)(0,1);
                    mu = irm_array.slice(j)(1,2);
                    
                    // compute and save the eigenvalues
                    real_eigenvals(0,j) = -beta;
                    real_eigenvals(1,j) = -mu;
                    
                    // compute the second and third eigenvectors
                    eigenvecs.slice(j)(0, 1) = beta / (beta - mu);

                    // compute the inverse eigenvector matrix
                    inversevecs.slice(j)(0,1) = -eigenvecs.slice(j)(0, 1);
                    inversevecs.slice(j)(0,2) = eigenvecs.slice(j)(0, 1) - 1;
          }
}

//' Update eigen values, vectors, and inverse matrices analytically for the SEIR 
//' model
//'
//' @param eigenvals matrix of eigenvalues with columns corresponding to irms
//' @param eigenvecs array of eigenvectors with third dimension corresponding to
//'   irms
//' @param inversevecs array of inverses of matrices of eigenvectors
//' @param irm_array array of rate matrices
//' @param complex_eigs logical vector indicating whether any eigen values are complex
//' @param first_calc boolean indicating whether this is the first time that the 
//'   eigen decomposition is being computed
//'
//' @return Updated eigenvalues, eigenvectors, and inverse matrices
//' @export
//' 
// [[Rcpp::export]]
void buildEigenArray_SEIR(arma::mat& real_eigenvals, arma::mat& imag_eigenvals, arma::cube& eigenvecs, arma::cube& inversevecs, arma::cube& irm_array, Rcpp::IntegerVector& n_real_eigs, bool initial_calc) {
          
          // Get dimensions of the irm array
          int num_irms = irm_array.n_slices;
          
          // fourth eigenvalue is always 0 and fourth eigenvector is always 1
          if(initial_calc) {
                    // zero out eigenvalues and eigenvectors
                    real_eigenvals.zeros();
                    imag_eigenvals.zeros();
                    eigenvecs.zeros();
                    inversevecs.zeros();
                    
                    // fill out eigenvector elements that are known
                    eigenvecs.each_slice() = arma::eye(4,4);      // set each each eigenvector matrix to diag(1,4)
                    eigenvecs.tube(0,3, arma::size(3, 1)).ones(); // 0 eigenvector is proportional to 1
                    
                    // fill out the elements of the inverse matrices that are known
                    inversevecs.each_slice() = arma::eye(4,4);    // set the diagonal of each inverse mtx to 1
                    inversevecs.tube(2,3).fill(-1); // set the (3,4) element to -1
                    
                    // all four eigenvectors are real
                    n_real_eigs.fill(4);
          }
          
          // initialize parameter values
          double beta  = 0;
          double gamma = 0;
          double mu    = 0;
          
          // Compute new eigen decompositions
          for(int j = 0; j < num_irms; ++j) {
                    
                    beta  = irm_array.slice(j)(0,1);
                    gamma = irm_array.slice(j)(1,2);
                    mu    = irm_array.slice(j)(2,3);
                    
                    // compute and save the eigenvalues
                    real_eigenvals(0,j) = -beta;
                    real_eigenvals(1,j) = -gamma;
                    real_eigenvals(2,j) = -mu;
                    
                    // compute the second eigenvector (for gamma)
                    eigenvecs.slice(j)(1, 1) = 1;
                    eigenvecs.slice(j)(0, 1) = beta / (beta - gamma);
                    
                    // compute the third eigenvector (for mu)
                    eigenvecs.slice(j)(2, 2) = 1;
                    eigenvecs.slice(j)(1, 2) = gamma / (gamma - mu);
                    eigenvecs.slice(j)(0, 2) = beta / (beta - mu) * eigenvecs.slice(j)(1, 2);
                    
                    // compute the inverse eigenvector matrix
                    inversevecs.slice(j)(0,1) = -eigenvecs.slice(j)(0,1);
                    inversevecs.slice(j)(1,2) = -eigenvecs.slice(j)(1,2);
                    inversevecs.slice(j)(0,2) = eigenvecs.slice(j)(0,1)*eigenvecs.slice(j)(1,2) - eigenvecs.slice(j)(0,2);
                    inversevecs.slice(j)(1,3) = eigenvecs.slice(j)(1,2) - 1;
                    inversevecs.slice(j)(0,3) = eigenvecs.slice(j)(0,1) - inversevecs.slice(j)(0,2) - 1;
          }
}

// //' Update eigen values, vectors, and inverse matrices analytically for the SIRS 
// //' model
// //'
// //' @param eigenvals matrix of eigenvalues with columns corresponding to irms
// //' @param eigenvecs array of eigenvectors with third dimension corresponding to
// //'   irms
// //' @param inversevecs array of inverses of matrices of eigenvectors
// //' @param irm_array array of rate matrices
// //' @param complex_eigs logical vector indicating whether any eigen values are complex
// //' @param first_calc boolean indicating whether this is the first time that the 
// //'   eigen decomposition is being computed
// //'
// //' @return Updated eigenvalues, eigenvectors, and inverse matrices
// //' @export
// //' 
// // [[Rcpp::export]]
// void buildEigenArray_SIRS(arma::mat& real_eigenvals, arma::mat& imag_eigenvals, arma::cube& eigenvecs, arma::cube& inversevecs, arma::cube& irm_array, Rcpp::IntegerVector& n_real_eigs, bool initial_calc) {
//           
//           // Get dimensions of the irm array
//           int num_rows = irm_array.n_rows;
//           int num_irms = irm_array.n_slices;
// 
//           // first eigenvalue is always 0 and first eigenvector is always 1
//           if(initial_calc) {
//                     real_eigenvals.zeros();
//                     imag_eigenvals.zeros();
//                     eigenvecs.zeros();
//                     inversevecs.zeros();
//                     eigenvecs.tube(0,0,arma::size(num_rows, 1)).ones(); // set the first eigenvector to 1
//                     eigenvecs.tube(0,1).ones();         // set the first element of the second eigenvector to 1
//           }
//           
//           // initialize parameter values
//           double beta = 0;
//           double mu = 0;
//           double gamma = 0;
//           double D = 0;
// 
//           // Compute new eigen decompositions
//           for(int j = 0; j < num_irms; ++j) {
//                     
//                     beta = irm_array.slice(j)(0,1);
//                     mu = irm_array.slice(j)(1,2);
//                     gamma = irm_array.slice(j)(2,0);
//                     D = pow((beta-mu),2.0) + gamma * (gamma - 2*(beta + mu));
// 
//                     if(D >= 0) {
//                               // set the number of real eigenvals
//                               n_real_eigs[j] = 3;
//                               
//                               // compute and save the eigenvalues
//                               real_eigenvals(1,j) = 0.5 * (-(beta + mu + gamma) + sqrt(D));
//                               real_eigenvals(2,j) = 0.5 * (-(beta + mu + gamma) - sqrt(D));
//                               imag_eigenvals(arma::span(1,2), j).zeros(); // set the imaginary parts of the second and third eigenvals to 0          
//                               
//                               // compute the second eigenvector 
//                               eigenvecs.slice(j)(2, 1) = gamma / (gamma + real_eigenvals(1,j));
//                               eigenvecs.slice(j)(1, 1) = mu / (mu + real_eigenvals(1,j)) * eigenvecs.slice(j)(2, 1);
//                               
//                               // compute the third eigenvector
//                               eigenvecs.slice(j)(0, 2) = 1;
//                               eigenvecs.slice(j)(2, 2) = gamma / (gamma + real_eigenvals(2,j));
//                               eigenvecs.slice(j)(1, 2) = mu / (mu + real_eigenvals(2,j)) * eigenvecs.slice(j)(2, 2);
//                               
//                               inversevecs.slice(j) = arma::inv(eigenvecs.slice(j));
//                               
//                     } else {
//                               // set the number of real eigenvals
//                               n_real_eigs[j] = 1;
//                               
//                               // compute the real and complex parts of the second eigenvalue
//                               double r_part  = -0.5 * (beta + mu + gamma);
//                               double i_part  = 0.5 * sqrt(-D);
//                               
//                               // compute the denominators: ((gamma + r_part)^2 + i_part^2), and ((mu + r_part)^2 + i_part^2)
//                               double i_part2   = -D/4;
//                               double den_gamma = pow(gamma + r_part, 2.0) + i_part2;
//                               double den_mu    = pow(mu + r_part, 2.0) + i_part2;
// 
//                               // compute and save the eigenvalues
//                               real_eigenvals(1,j) = r_part;
//                               real_eigenvals(2,j) = r_part;
//                               imag_eigenvals(1,j) = i_part;
//                               imag_eigenvals(2,j) = -i_part;
//                               
//                               // The first element of the complex part of the second eigenvector is zero
//                               eigenvecs.slice(j)(0, 2) = 0;
//                               
//                               // compute the real and complex parts of the third element of the second eigenvector
//                               eigenvecs.slice(j)(2, 1) = gamma * (gamma + r_part) / den_gamma;
//                               eigenvecs.slice(j)(2, 2) = -gamma * i_part / den_gamma;
//                               
//                               // compute the real and complex parts of the second element of the second eigenvector
//                               eigenvecs.slice(j)(1, 1) = mu * ((mu+r_part)*eigenvecs.slice(j)(2,1) + eigenvecs.slice(j)(2,2)*i_part) / den_mu;
//                               eigenvecs.slice(j)(1, 2) = mu * ((mu+r_part)*eigenvecs.slice(j)(2,2) - eigenvecs.slice(j)(2,1)*i_part) / den_mu;
//                               
//                               // compute the inverse eigenvector matrix
//                               inversevecs.slice(j) = arma::inv(eigenvecs.slice(j));
//                     }
//           }
// }