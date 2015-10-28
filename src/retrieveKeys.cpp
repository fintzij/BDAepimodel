// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Get the irm keys for the compartment counts in a population level
//' bookkeeping matrix
//'
//' @param inds vector of indices, generally 1:epimodel$ind_final_config
//' @param irm_lookup lookup matrix relating configurations to irm keys
//' @param pop_mat population bookkeeping matrix in epimodel list
//' @param index_state_number vector indicating which states are index states
//'
//' @return Vector of keys to index into an IRM array.
// [[Rcpp::export]]
arma::vec retrieveKeys(Rcpp::IntegerVector inds, const arma::mat& irm_lookup, const arma::mat& pop_mat, const arma::vec& index_state_num) {

        // ensure that inds starts at 0 so subtract 1
        inds = inds - 1;

        // Get relevant dimensions
        int n_configs = inds.size(); // length of 1:ind_final_config
        int n_keys = irm_lookup.n_rows; // number of keys

        // Initialize output vector
        arma::vec keys(n_configs, arma::fill::zeros);

        // Initialize vector for tracking which configurations have been matched
        arma::vec config_inds(n_configs, arma::fill::zeros);

        // Get columns in pop_mat and irm_lookup to compare
        arma::vec pop_index_num(index_state_num + 2); // index state configurations start after three columns
        arma::uvec pop_index_inds = arma::conv_to<arma::uvec>::from(pop_index_num);

        arma::uvec lookup_index_inds(2);
        lookup_index_inds[0] = 1; // index state configurations start after one column
        lookup_index_inds[1] = index_state_num.size();

        // matrices corresponding to index states
        arma::mat p_config = pop_mat.cols(pop_index_inds); // only the columns from the pop_mat relating to the index states
        arma::mat l_config = irm_lookup.cols(arma::span(lookup_index_inds[0], lookup_index_inds[1]));

        // Fill out the keys
        // outer loop over unique configurations
        for(int k = 0; k < n_keys; ++k) {

                // Get config in the lookup matrix
                arma::rowvec l_conf = l_config.row(k);

                // inner loop over unmatched indices in pop_mat
                for(int j = 0; j < n_configs; ++j) {

                        // If the configuration has been matched
                        if(config_inds[j] != 1) {

                                arma::rowvec p_conf = p_config.row(j);

                                if(all(p_conf == l_conf)) {

                                        keys[j] = irm_lookup(k, 0); // assign the key
                                        config_inds[k] = 1; // indicate that the configuration has been matched
                                }
                        }
                }

                if(all(config_inds) == 1) break;
        }

        return keys;
}