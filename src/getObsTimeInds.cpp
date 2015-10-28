#include <Rcpp.h>
using namespace Rcpp;

//' Get the indices of observation times in the population level bookkeeping mtx.
//'
//' @param pop_mat population level bookkeeping matrix
//' @param obstimes vector of observation times
//'
//' @return indices of observation times
// [[Rcpp::export]]
Rcpp::IntegerVector getObsTimeInds(const Rcpp::NumericMatrix& pop_mat, const Rcpp::NumericVector& obstimes) {

        int pop_nrow = pop_mat.nrow();
        int n_obs = obstimes.size();

        Rcpp::IntegerVector inds(n_obs);
        int k = 1;

        inds[0] = 1;

        for(int j=1; j < pop_nrow; ++j) {

                if(pop_mat(j, 0) == obstimes[k]) {
                        inds[k] = j+1;
                        k += 1;
                }
                if(k == n_obs) break;
        }

        return inds;
}

