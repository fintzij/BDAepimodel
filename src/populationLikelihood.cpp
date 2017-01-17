// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

//' Get the irm keys for the compartment counts in a population level
//' bookkeeping matrix
//'
//' @param pop_mat population level bookkeeping matrix
//' @param irm array of rate matrices
//' @param initdist vector of initial state probabilities
//' @param initdist_param_inds vector of indices for admissible initial states
//' @param flow_inds matrix of indices for locations in IRM for each event
//' @param keys vector of irm keys
//' @param inds indices relating to event times, along with t0 and tmax
//' @param loglik boolean indicating whether to return the likelihood or
//'     log-likelihood
//'
//' @return population level likelihood or log-likelihood
// [[Rcpp::export]]
double populationLikelihood(const arma::mat& pop_mat, const arma::cube& irm, const arma::vec& initdist, const arma::uvec& initdist_param_inds, const arma::umat& flow_inds, const arma::uvec& keys, const arma::uvec& inds, bool loglik) {
          
          // get number of endpoints and number of intervals
          arma::uword n_intervals = inds.n_elem - 1;
          
          // ensure that indices (excluding obs_time_inds) start at 0
          arma::uvec pop_inds = inds - 1;
          arma::uvec init_inds = initdist_param_inds - 1;
          arma::uvec non_zero_inits = arma::find(initdist(init_inds) != 0);
          
          arma::umat event_inds = flow_inds - 1;
          arma::uword next_event(0);
          
          // get the relevant vectors
          arma::vec times_all = pop_mat.col(0);
          arma::vec event_seq_all = pop_mat.col(2);
          arma::mat counts_all = pop_mat(arma::span::all, span(3,pop_mat.n_cols - 1));

          // subset the vectors and count matrix so they only include event
          // times and the endpoints of the observation period
          arma::vec timediffs = diff(times_all.elem(pop_inds));
          arma::vec event_seq = event_seq_all.elem(pop_inds) - 1;
          arma::mat counts = counts_all.rows(pop_inds);
          arma::uvec irm_keys = keys.elem(pop_inds) - 1;

          // initialize the likelihood
          arma::mat init_counts = counts.row(0);

          // get indices of any initial state probabilities that equal zero
          double pop_lik = sum(init_counts(non_zero_inits) % log(initdist(non_zero_inits)));
          
          // compute the likelihood of the path
          for(uword j = 0; j < n_intervals - 1; ++j) {
                    next_event = event_seq[j+1];
                    pop_lik += log(irm(event_inds(next_event, 0), event_inds(next_event, 1), irm_keys[j])) + arma::accu(irm.slice(irm_keys[j]).diag() % counts.row(j).t()) * timediffs[j];
          }
          
          // Tail probability
          pop_lik += arma::accu(irm.slice(irm_keys[n_intervals]).diag() % counts.row(n_intervals).t()) * timediffs[n_intervals - 1];
          
          if(loglik) {
                    return pop_lik;
          } else {
                    return exp(pop_lik);
          }
}