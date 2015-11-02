// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Insert subject transitions into the population-level and subject level 
//' bookkeeping matrices.
//'
//' @param path matrix of transitions to be inserted
//' @param subject subject ID
//' @param pop_mat population bookkeeping matrix in epimodel list
//' @param subj_path vector containing the subject path
//' @param subj_row_ind index stating where to begin inserting the transitions
//'
//' @return updated bookeeping objects
// [[Rcpp::export]]
void insertPath(const Rcpp::NumericMatrix& path, const int subject, Rcpp::NumericMatrix& pop_mat, Rcpp::IntegerVector& subj_path, int& ind) {
          
          // Get path dimensions
          int n_jumps = path.nrow();

          // insert path
          for(int j = 0; j < n_jumps; ++j) {

                    pop_mat(ind + j, 0) = path(j, 0); // insert time
                    pop_mat(ind + j, 1) = subject; // insert subject ID
                    pop_mat(ind + j, 2) = path(j, 1); // insert event ID
                    
                    subj_path(ind + j) = int(path(j, 2)); // insert state into config mat
          }
}