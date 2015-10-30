// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Get the irm keys for the compartment counts in a population level
//' bookkeeping matrix
//'
//' @param subject
//' @param pop_mat population level bookkeeping matrix
//' @param config_mat subject_level bookkeeping matrix for configurations
//' @param irm_array array of rate matrices
//' @param keys vector of irm keys
//' @param inds indices relating to event times, along with t0 and tmax
//' @param loglik boolean indicating whether to return the likelihood or
//'     log-likelihood
//'
//' @return subject level likelihood or log-likelihood
// [[Rcpp::export]]
double subjectLikelihood(const int subject, const arma::mat& pop_mat, const arma::mat& config_mat, Rcpp::NumericVector& irm_array, const Rcpp::NumericVector& initdist, const arma::uvec& keys, const arma::uvec& inds, bool loglik) {

        int n_endpoints = inds.n_elem;
        arma::uword n_intervals = n_endpoints - 1;

        arma::uvec subj_inds = inds - 1;

        // Instatiate array pointers
        Rcpp::IntegerVector irmDims = irm_array.attr("dim");
        arma::cube irm(irm_array.begin(), irmDims[0], irmDims[1], irmDims[2], false);

        // get the relevant vectors
        arma::vec times_all = pop_mat.col(0);
        arma::vec subj_path_all = config_mat.col(subject - 1);

        // subset the vectors so they only include event times and the endpoints of the observation period
        arma::vec timediffs = diff(times_all.elem(subj_inds));
        arma::vec subj_path = subj_path_all.elem(subj_inds) - 1;
        arma::uvec irm_keys = keys.elem(subj_inds) - 1;

        // initialize the likelihood
        double subj_lik = log(initdist[subj_path[0]]);

        // initialize the state and dt variables
        int next_state = subj_path[0];
        int cur_state = subj_path[0];

        // compute the likelihood of the path
        for(uword j = 0; j < n_intervals; ++j) {

                // increment the states
                cur_state = next_state;
                next_state = subj_path[j+1];

                if(next_state == cur_state) {

                        subj_lik += irm(cur_state, cur_state, irm_keys[j]) * timediffs[j];


                } else {

                        subj_lik += log(irm(cur_state, next_state, irm_keys[j])) + irm(cur_state, cur_state, irm_keys[j]) * timediffs[j];

                }

        }

        if(loglik) {
                return subj_lik;
        } else {
                return exp(subj_lik);
        }
}