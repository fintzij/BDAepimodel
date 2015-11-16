// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

using namespace arma;
using namespace Rcpp;

//' Rcpp code to simulate a general stochastic epidemic and binomial counts. 
//'
//' @param popsize population size
//' @param obstimes vector of observation times
//' @param params vector of parameters: beta, mu
//' @param init_config initial configuration of compartment counts
//'
//' @return updated bookeeping objects
// [[Rcpp::export]]
Rcpp::NumericMatrix simulateSIR(Rcpp::NumericVector obstimes, Rcpp::NumericVector params, Rcpp::IntegerVector init_config) {
          
          // initialize population bookkeeping matrix
          int popsize = sum(init_config);
          Rcpp::NumericMatrix pop_mat(2*popsize, 6);
          pop_mat(0, 3) = init_config[0];
          pop_mat(0, 4) = init_config[1];
          pop_mat(0, 5) = init_config[2];
          
          // get parameters
          double beta = params[0]; 
          double mu = params[1];
          
          // set initial values
          int nobs = obstimes.size();
          double t = obstimes[0]; // current time
          Rcpp::NumericVector dt(1); // time increment
          double tmax = obstimes[nobs - 1]; // end of observation period

          double St = init_config[0]; // number of susceptibles
          double It = init_config[1]; // number of infecteds
          double Rt = init_config[2]; // number of recovereds
          
          // Vector of lumped rates
          Rcpp::NumericVector rates(2);
          rates[0] = beta*St*It;
          rates[1] = mu*It; 
          
          Rcpp::IntegerVector events = Rcpp::seq_len(2); // vector of event codes
          Rcpp::IntegerVector next_event(1); // next event
          Rcpp::IntegerVector subj_ind(1); // subject index
          int subject(1); // subject ID
          
          Rcpp::IntegerVector susc_IDs = Rcpp::seq_len(St); // IDs for susceptibles
          Rcpp::IntegerVector infec_IDs = St + Rcpp::seq_len(It); // IDs for infecteds
          Rcpp::IntegerVector susc_Inds = Rcpp::seq_len(St) - 1; // inds for susceptibles
          Rcpp::IntegerVector infec_Inds = Rcpp::seq_len(It) - 1; // inds for infecteds
          std::vector<int> susceptibles = Rcpp::as<std::vector<int> >(susc_IDs); // convert to std vector
          std::vector<int> infecteds = Rcpp::as<std::vector<int> >(infec_IDs);  // convert to std vector
          std::vector<int> susc_inds = Rcpp::as<std::vector<int> >(susc_Inds); // convert to std vector
          std::vector<int> infec_inds = Rcpp::as<std::vector<int> >(infec_Inds);  // convert to std vector
          
          // set keep_going and the row index
          bool keep_going = true;
          int ind = 0;

          // start simulating
          while(keep_going) {
                    
                    // sample the next time
                    dt = Rcpp::rexp(1, sum(rates));
                    t += dt[0];
                    
                    if(t > tmax) {
                              keep_going = false; // stop simulating
                              
                    } else {
                              // sample the next event
                              next_event = Rcpp::RcppArmadillo::sample(events, 1, false, rates);
                              
                              if(next_event[0] == 1) { // infection occurs
                                        
                                        // choose subject
                                        subject = susceptibles.back();
                                        
                                        // move subject
                                        infecteds.push_back(subject); // move to infecteds
                                        susceptibles.pop_back(); // remove from susc IDs
                                        
                                        susc_inds.pop_back(); // shrink susceptible indices vector
                                        infec_inds.push_back(infec_inds.size()); // expand infected indices
                                        
                                        // update compartment counts
                                        St -= 1;
                                        It += 1;
                                        
                                        // update rates
                                        rates[0] = beta * St * It;
                                        rates[1] = mu * It; 
                                        
                                        // update bookkeeping matrix
                                        pop_mat(ind, 0) = t; // new time
                                        pop_mat(ind, 1) = subject; // subject ID
                                        pop_mat(ind, 2) = 1; // event code
                                        pop_mat(ind, 3) = St; // number of susceptibles
                                        pop_mat(ind, 4) = It; // number of infecteds
                                        pop_mat(ind, 5) = Rt; // number of recovereds
                                        
                                        // increment row index
                                        ind += 1;
                                        
                                        // check if any more transitions are possible (i.e. all rates equal to zero).
                                        // if not, set keep_going to false.
                                        if((rates[0] == 0) & (rates[1] == 0)) {
                                                  keep_going = false;
                                        }
                                        
                              } else if(next_event[0] == 2) { // recovery occurs
                                        
                                        // choose subject
                                        subj_ind = Rcpp::RcppArmadillo::sample(infec_inds, 1, false, Rcpp::NumericVector::create());
                                        subject = infecteds[subj_ind[0]];
                                        
                                        // move subject
                                        infecteds.erase(infecteds.begin() + subj_ind[0]); // remove from susc IDs
                                        infec_inds.pop_back(); // shrink infected indices
                                        
                                        // update compartment counts
                                        It -= 1;
                                        Rt += 1;
                                        
                                        // update rates
                                        rates[0] = beta * St * It;
                                        rates[1] = mu * It; 
                                        
                                        // update bookkeeping matrix
                                        pop_mat(ind, 0) = t; // new time
                                        pop_mat(ind, 1) = subject; // subject ID
                                        pop_mat(ind, 2) = 2; // event code
                                        pop_mat(ind, 3) = St; // number of susceptibles
                                        pop_mat(ind, 4) = It; // number of infecteds
                                        pop_mat(ind, 5) = Rt; // number of recovereds
                                        
                                        // increment row index
                                        ind += 1;
                                        
                                        // check if any more transitions are possible (i.e. all rates equal to zero).
                                        // if not, set keep_going to false.
                                        if((rates[0] == 0) & (rates[1] == 0)) {
                                                  keep_going = false;
                                        }
                                        
                              }
                    }
          }
          
          return pop_mat;
}