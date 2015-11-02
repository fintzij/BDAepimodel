// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Insert subject transitions into the population-level and subject level 
//' bookkeeping matrices.
//'
//' @param subj_path subject path vector
//' @param subject ID for subject
//' @param pop_mat population level bookkeeping matrix
//' @param init_config initial configuration vector
//' @param ind_final_config index for final configuration
//' @param flow_inds matrix of flow indices
//'
//' @return extended path vector
// [[Rcpp::export]]
void retrieveSubjPath(Rcpp::IntegerVector& subj_path, const int subject, const Rcpp::NumericMatrix& pop_mat, const Rcpp::IntegerVector& init_config, const int ind_final_config, const Rcpp::IntegerMatrix flow_inds) {
          
          // Get object dimensions
          Rcpp::IntegerVector popDims = pop_mat.attr("dim");

          // initialize state
          int cur_state(1);
          int event(1);
          
          // set initial state
          subj_path[0] = cur_state = init_config[subject - 1];
          
          // record the states
          for(int j = 1; j < ind_final_config; ++j) {
                    
                    // if the next event did not belong to the subject, record the state
                    if(pop_mat(j, 1) != subject) {
                              
                              subj_path[j] = cur_state; // record the state

                    } else {
                              event = pop_mat(j, 2) - 1; // get the row in the flow_inds matrix
                              subj_path[j] = cur_state = flow_inds(event, 1); // get the next state
                    }
          }
}