// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

//' Update eigen values, vectors, and inverse matrices for irms
//'
//' @param tpms array of transition probability matrices to be modified
//' @param pop_mat population level bookkeeping matrix
//' @param eigen_vals matrix of eigen values of irms
//' @param eigen_vecs array of eigen vectors of rate matrices
//' @param inverse_vecs array of inverse matrices of eigen vectors
//' @param irm_keys vector of irm array keys
//'
//' @return Updated eigenvalues, eigenvectors, and inverse matrices
// [[Rcpp::export]]
Rcpp::IntegerVector sampleEventSubseq(Rcpp::IntegerVector& path, Rcpp::NumericVector& tpms, Rcpp::NumericVector& tpm_prods, const int init_ind, const int final_ind) {

          // set the initial and final inds
          int ind_start = init_ind - 1;
          int ind_end = final_ind -1;

          // set the initial and final states
          int state_cur = path[ind_start];
          int state_end = path[ind_end];
          Rcpp::IntegerVector next_state(1);

          // get tpm dimensions
          Rcpp::IntegerVector tpmDims = tpms.attr("dim");

          // set tpm pointers and path pointer
          arma::cube tpm_arr(tpms.begin(), tpmDims[0], tpmDims[1], tpmDims[2], false);
          arma::cube prod_arr(tpm_prods.begin(), tpmDims[0], tpmDims[1], tpmDims[2], false);

          // create vector for state probabilities
          Rcpp::NumericVector state_probs(tpmDims[0]);
          Rcpp::IntegerVector states = Rcpp::seq_len(tpmDims[0]);

          // Sample the status at observation times
          for(int j = (ind_start+1); j < ind_end; ++j) {

                    state_probs = (tpm_arr.slice(j - 1).row(state_cur - 1).t() % prod_arr.slice(j).col(state_end - 1)) / prod_arr.slice(j - 1).at(state_cur - 1, state_end - 1);

                    next_state = Rcpp::RcppArmadillo::sample(states, 1, false, state_probs);

                    path[j] = state_cur = int(next_state[0]);
          }

          return path;
}

//' simulateSIR function, included in this file because of difficulty including sample.h in two separate files.
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
          Rcpp::NumericVector probs(2);
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
                              probs = rates / sum(rates);
                              next_event = Rcpp::RcppArmadillo::sample(events, 1, false, probs);

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
// 
// //' simulateSEIR function, included in this file because of difficulty including sample.h in two separate files.
// //' Rcpp code to simulate a general stochastic epidemic and binomial counts.
// //'
// //' @param popsize population size
// //' @param obstimes vector of observation times
// //' @param params vector of parameters: beta, mu
// //' @param init_config initial configuration of compartment counts
// //'
// //' @return updated bookeeping objects
// // [[Rcpp::export]]
// Rcpp::NumericMatrix simulateSEIR(Rcpp::NumericVector obstimes, Rcpp::NumericVector params, Rcpp::IntegerVector init_config) {
//           
//           // initialize population bookkeeping matrix
//           int popsize = sum(init_config);
//           Rcpp::NumericMatrix pop_mat(3*popsize, 7);
//           pop_mat(0, 3) = init_config[0];
//           pop_mat(0, 4) = init_config[1];
//           pop_mat(0, 5) = init_config[2];
//           pop_mat(0, 6) = init_config[3];
//           
//           // get parameters
//           double beta  = params[0];
//           double gamma = params[1];
//           double mu    = params[2];
//           
//           // set initial values 
//           int nobs = obstimes.size();
//           double t = obstimes[0];           // current time
//           Rcpp::NumericVector dt(1);        // time increment
//           double tmax = obstimes[nobs - 1]; // end of observation period
//           
//           double St = init_config[0]; // number of susceptibles
//           double Et = init_config[1]; // number of exposed
//           double It = init_config[2]; // number of infectious
//           double Rt = init_config[3]; // number of recovereds
//           
//           // Vector of lumped rates
//           Rcpp::NumericVector rates(3);
//           Rcpp::NumericVector probs(3);
//           rates[0] = beta*St*It;
//           rates[1] = gamma*Et;
//           rates[2] = mu*It;
//           
//           Rcpp::IntegerVector events = Rcpp::seq_len(3); // vector of event codes
//           Rcpp::IntegerVector next_event(1);             // next event
//           Rcpp::IntegerVector subj_ind(1);               // subject index
//           int subject(1);                                // subject ID
//           
//           Rcpp::IntegerVector susc_IDs   = Rcpp::seq_len(St);            // IDs for susceptibles
//           Rcpp::IntegerVector exp_IDs    = St + Rcpp::seq_len(Et);       // IDs for exposed
//           Rcpp::IntegerVector infec_IDs  = St + Et + Rcpp::seq_len(It);  // IDs for infecteds
//           Rcpp::IntegerVector susc_Inds  = Rcpp::seq_len(St) - 1;        // inds for susceptibles
//           Rcpp::IntegerVector exp_Inds   = Rcpp::seq_len(Et) - 1;        // inds for exposed
//           Rcpp::IntegerVector infec_Inds = Rcpp::seq_len(It) - 1;        // inds for infecteds
//           std::vector<int> susceptibles  = Rcpp::as<std::vector<int> >(susc_IDs);     // convert to std vector
//           std::vector<int> exposed       = Rcpp::as<std::vector<int> >(exp_IDs);      // convert to std vector
//           std::vector<int> infecteds     = Rcpp::as<std::vector<int> >(infec_IDs);    // convert to std vector
//           std::vector<int> susc_inds     = Rcpp::as<std::vector<int> >(susc_Inds);    // convert to std vector
//           std::vector<int> exp_inds      = Rcpp::as<std::vector<int> >(exp_Inds);     // convert to std vector
//           std::vector<int> infec_inds    = Rcpp::as<std::vector<int> >(infec_Inds);   // convert to std vector
//           
//           // set keep_going and the row index
//           bool keep_going = true;
//           int ind = 0;
//           
//           // start simulating
//           while(keep_going) {
//                     
//                     // sample the next time
//                     dt = Rcpp::rexp(1, sum(rates));
//                     t += dt[0];
//                     
//                     if(t > tmax) {
//                               keep_going = false; // stop simulating
//                               
//                     } else {
//                               // sample the next event
//                               probs = rates / sum(rates);
//                               next_event = Rcpp::RcppArmadillo::sample(events, 1, false, probs);
//                               
//                               if(next_event[0] == 1) { // exposure occurs
//                                         
//                                         // choose subject
//                                         subject = susceptibles.back();
//                                         
//                                         // move subject
//                                         exposed.push_back(subject); // move to infecteds
//                                         susceptibles.pop_back();    // remove from susc IDs
//                                         
//                                         susc_inds.pop_back();                // shrink susceptible indices vector
//                                         exp_inds.push_back(exp_inds.size()); // expand infected indices
//                                         
//                                         // update compartment counts
//                                         St -= 1;
//                                         Et += 1;
//                                         
//                                         // update rates
//                                         rates[0] = beta * St * It;
//                                         rates[1] = gamma * Et;
//                                         rates[2] = mu * It;
// 
//                                         // update bookkeeping matrix
//                                         pop_mat(ind, 0) = t; // new time
//                                         pop_mat(ind, 1) = subject; // subject ID
//                                         pop_mat(ind, 2) = 1; // event code
//                                         pop_mat(ind, 3) = St; // number of susceptibles
//                                         pop_mat(ind, 4) = Et; // number of exposed
//                                         pop_mat(ind, 5) = It; // number of infecteds
//                                         pop_mat(ind, 6) = Rt; // number of recovereds
//                                         
//                                         // increment row index
//                                         ind += 1;
//                                         
//                                         // check if any more transitions are possible (i.e. all rates equal to zero).
//                                         // if not, set keep_going to false.
//                                         if((rates[0] == 0) & (rates[1] == 0) & (rates[2] == 0)) {
//                                                   keep_going = false;
//                                         }
//                                         
//                               } else if(next_event[0] == 2) { // subject becomes infectious
//                                         
//                                         // choose subject
//                                         subj_ind = Rcpp::RcppArmadillo::sample(exp_inds, 1, false, Rcpp::NumericVector::create());
//                                         subject = exposed[subj_ind[0]];
//                                         
//                                         // move subject
//                                         exposed.erase(exposed.begin() + subj_ind[0]); // remove from exposed IDs
//                                         infecteds.push_back(subject);                 // move to infecteds
//                                         
//                                         exp_inds.pop_back();                          // shrink infected indices
//                                         infec_inds.push_back(infec_inds.size());      // expand infected indices
//                                         
//                                         // update compartment counts
//                                         Et -= 1;
//                                         It += 1;
//                                         
//                                         // update rates
//                                         rates[0] = beta * St * It;
//                                         rates[1] = gamma * Et;
//                                         rates[2] = mu * It;
//                                         
//                                         // update bookkeeping matrix
//                                         pop_mat(ind, 0) = t; // new time
//                                         pop_mat(ind, 1) = subject; // subject ID
//                                         pop_mat(ind, 2) = 2; // event code
//                                         pop_mat(ind, 3) = St; // number of susceptibles
//                                         pop_mat(ind, 4) = Et;
//                                         pop_mat(ind, 5) = It; // number of infecteds
//                                         pop_mat(ind, 6) = Rt; // number of recovereds
//                                         
//                                         // increment row index
//                                         ind += 1;
//                                         
//                                         // check if any more transitions are possible (i.e. all rates equal to zero).
//                                         // if not, set keep_going to false.
//                                         if((rates[0] == 0) & (rates[1] == 0) & (rates[2] == 0)) {
//                                                   keep_going = false;
//                                         }
//                                         
//                               } else if(next_event[0] == 3) { // subject becomes infectious
//                                         
//                                         // choose subject
//                                         subj_ind = Rcpp::RcppArmadillo::sample(infec_inds, 1, false, Rcpp::NumericVector::create());
//                                         subject = infecteds[subj_ind[0]];
//                                         
//                                         // move subject
//                                         infecteds.erase(infecteds.begin() + subj_ind[0]); // remove from infected IDs
//                                         infec_inds.pop_back(); // shrink infected indices
//                                         
//                                         // update compartment counts
//                                         It -= 1;
//                                         Rt += 1;
//                                         
//                                         // update rates
//                                         rates[0] = beta * St * It;
//                                         rates[1] = gamma * Et;
//                                         rates[2] = mu * It;
//                                         
//                                         // update bookkeeping matrix
//                                         pop_mat(ind, 0) = t; // new time
//                                         pop_mat(ind, 1) = subject; // subject ID
//                                         pop_mat(ind, 2) = 3; // event code
//                                         pop_mat(ind, 3) = St; // number of susceptibles
//                                         pop_mat(ind, 4) = Et;
//                                         pop_mat(ind, 5) = It; // number of infecteds
//                                         pop_mat(ind, 6) = Rt; // number of recovereds
//                                         
//                                         // increment row index
//                                         ind += 1;
//                                         
//                                         // check if any more transitions are possible (i.e. all rates equal to zero).
//                                         // if not, set keep_going to false.
//                                         if((rates[0] == 0) & (rates[1] == 0) & (rates[2] == 0)) {
//                                                   keep_going = false;
//                                         }
//                               }
//                     }
//           }
//           
//           return pop_mat;
// }
// 
// //' simulateSIIR function, included in this file because of difficulty including sample.h in two separate files.
// //' Rcpp code to simulate a general stochastic epidemic and binomial counts.
// //'
// //' @param popsize population size
// //' @param obstimes vector of observation times
// //' @param params vector of parameters: beta, mu
// //' @param init_config initial configuration of compartment counts
// //'
// //' @return updated bookeeping objects
// // [[Rcpp::export]]
// Rcpp::NumericMatrix simulateSIIR(Rcpp::NumericVector obstimes, Rcpp::NumericVector params, Rcpp::IntegerVector init_config) {
//           
//           // initialize population bookkeeping matrix
//           int popsize = sum(init_config);
//           Rcpp::NumericMatrix pop_mat(3*popsize, 7);
//           pop_mat(0, 3) = init_config[0];
//           pop_mat(0, 4) = init_config[1];
//           pop_mat(0, 5) = init_config[2];
//           pop_mat(0, 6) = init_config[3];
//           
//           // get parameters
//           double beta  = params[0];
//           double gamma = params[1];
//           double mu    = params[2];
//           
//           // set initial values 
//           int nobs = obstimes.size();
//           double t = obstimes[0];           // current time
//           Rcpp::NumericVector dt(1);        // time increment
//           double tmax = obstimes[nobs - 1]; // end of observation period
//           
//           double St = init_config[0];  // number of susceptibles
//           double I1t = init_config[1]; // number of infecteds, comp 1
//           double I2t = init_config[2]; // number of infectious, comp 2
//           double Rt = init_config[3];  // number of recovereds
//           
//           // Vector of lumped rates
//           Rcpp::NumericVector rates(3);
//           Rcpp::NumericVector probs(3);
//           rates[0] = beta*St*(I1t + I2t);
//           rates[1] = gamma*I1t;
//           rates[2] = mu*I2t;
//           
//           Rcpp::IntegerVector events = Rcpp::seq_len(3); // vector of event codes
//           Rcpp::IntegerVector next_event(1);             // next event
//           Rcpp::IntegerVector subj_ind(1);               // subject index
//           int subject(1);                                // subject ID
//           
//           Rcpp::IntegerVector susc_IDs    = Rcpp::seq_len(St);             // IDs for susceptibles
//           Rcpp::IntegerVector infec1_IDs  = St + Rcpp::seq_len(I1t);       // IDs for infecteds comp 1
//           Rcpp::IntegerVector infec2_IDs  = St + I1t + Rcpp::seq_len(I2t); // IDs for infecteds comp 2
//           Rcpp::IntegerVector susc_Inds   = Rcpp::seq_len(St) - 1;         // inds for susceptibles
//           Rcpp::IntegerVector infec1_Inds = Rcpp::seq_len(I1t) - 1;        // inds for exposed
//           Rcpp::IntegerVector infec2_Inds = Rcpp::seq_len(I2t) - 1;        // inds for infecteds
//           std::vector<int> susceptibles   = Rcpp::as<std::vector<int> >(susc_IDs);    // convert to std vector
//           std::vector<int> infecteds1     = Rcpp::as<std::vector<int> >(infec1_IDs);  // convert to std vector
//           std::vector<int> infecteds2     = Rcpp::as<std::vector<int> >(infec2_IDs);  // convert to std vector
//           std::vector<int> susc_inds      = Rcpp::as<std::vector<int> >(susc_Inds);   // convert to std vector
//           std::vector<int> infec1_inds    = Rcpp::as<std::vector<int> >(infec1_Inds); // convert to std vector
//           std::vector<int> infec2_inds    = Rcpp::as<std::vector<int> >(infec2_Inds); // convert to std vector
//           
//           // set keep_going and the row index
//           bool keep_going = true;
//           int ind = 0;
//           
//           // start simulating
//           while(keep_going) {
//                     
//                     // sample the next time
//                     dt = Rcpp::rexp(1, sum(rates));
//                     t += dt[0];
//                     
//                     if(t > tmax) {
//                               keep_going = false; // stop simulating
//                               
//                     } else {
//                               // sample the next event
//                               probs = rates / sum(rates);
//                               next_event = Rcpp::RcppArmadillo::sample(events, 1, false, probs);
//                               
//                               if(next_event[0] == 1) { // exposure occurs
//                                         
//                                         // choose subject
//                                         subject = susceptibles.back();
//                                         
//                                         // move subject
//                                         infecteds1.push_back(subject); // move to infecteds
//                                         susceptibles.pop_back();       // remove from susc IDs
//                                         
//                                         susc_inds.pop_back(); // shrink susceptible indices vector
//                                         infec1_inds.push_back(infec1_inds.size()); // expand infected indices
//                                         
//                                         // update compartment counts
//                                         St -= 1;
//                                         I1t += 1;
//                                         
//                                         // update rates
//                                         rates[0] = beta * St * (I1t + I2t);
//                                         rates[1] = gamma * I1t;
//                                         
//                                         // update bookkeeping matrix
//                                         pop_mat(ind, 0) = t; // new time
//                                         pop_mat(ind, 1) = subject; // subject ID
//                                         pop_mat(ind, 2) = 1; // event code
//                                         pop_mat(ind, 3) = St; // number of susceptibles
//                                         pop_mat(ind, 4) = I1t; // number of exposed
//                                         pop_mat(ind, 5) = I2t; // number of infecteds
//                                         pop_mat(ind, 6) = Rt; // number of recovereds
//                                         
//                                         // increment row index
//                                         ind += 1;
//                                         
//                                         // check if any more transitions are possible (i.e. all rates equal to zero).
//                                         // if not, set keep_going to false.
//                                         if((rates[0] == 0) & (rates[1] == 0) & (rates[2] == 0)) {
//                                                   keep_going = false;
//                                         }
//                                         
//                               } else if(next_event[0] == 2) { // subject becomes infectious
//                                         
//                                         // choose subject
//                                         subj_ind = Rcpp::RcppArmadillo::sample(infec1_inds, 1, false, Rcpp::NumericVector::create());
//                                         subject = infecteds1[subj_ind[0]];
//                                         
//                                         // move subject
//                                         infecteds1.erase(infecteds1.begin() + subj_ind[0]); // remove from infecteds 1
//                                         infecteds2.push_back(subject);                      // move to infecteds 2
//                                         
//                                         infec1_inds.pop_back();                    // shrink infected 1 indices
//                                         infec2_inds.push_back(infec2_inds.size()); // expand infected 2 indices
//                                         
//                                         // update compartment counts
//                                         I1t -= 1;
//                                         I2t += 1;
//                                         
//                                         // update rates
//                                         rates[0] = beta * St * (I1t + I2t);
//                                         rates[1] = gamma * I1t;
//                                         rates[2] = mu * I2t;
//                                         
//                                         // update bookkeeping matrix
//                                         pop_mat(ind, 0) = t; // new time
//                                         pop_mat(ind, 1) = subject; // subject ID
//                                         pop_mat(ind, 2) = 2;   // event code
//                                         pop_mat(ind, 3) = St;  // number of susceptibles
//                                         pop_mat(ind, 4) = I1t; // number of infecteds 1
//                                         pop_mat(ind, 5) = I2t; // number of infecteds 2
//                                         pop_mat(ind, 6) = Rt;  // number of recovereds
//                                         
//                                         // increment row index
//                                         ind += 1;
//                                         
//                                         // check if any more transitions are possible (i.e. all rates equal to zero).
//                                         // if not, set keep_going to false.
//                                         if((rates[0] == 0) & (rates[1] == 0) & (rates[2] == 0)) {
//                                                   keep_going = false;
//                                         }
//                                         
//                               } else if(next_event[0] == 3) { // subject becomes infectious
//                                         
//                                         // choose subject
//                                         subj_ind = Rcpp::RcppArmadillo::sample(infec2_inds, 1, false, Rcpp::NumericVector::create());
//                                         subject = infecteds2[subj_ind[0]];
//                                         
//                                         // move subject
//                                         infecteds2.erase(infecteds2.begin() + subj_ind[0]); // remove from infected IDs
//                                         infec2_inds.pop_back(); // shrink infected indices
//                                         
//                                         // update compartment counts
//                                         I2t -= 1;
//                                         Rt  += 1;
//                                         
//                                         // update rates
//                                         rates[0] = beta * St * (I1t + I2t);
//                                         rates[2] = mu * I2t;
//                                         
//                                         // update bookkeeping matrix
//                                         pop_mat(ind, 0) = t;   // new time
//                                         pop_mat(ind, 1) = subject; // subject ID
//                                         pop_mat(ind, 2) = 3;   // event code
//                                         pop_mat(ind, 3) = St;  // number of susceptibles
//                                         pop_mat(ind, 4) = I1t; // number of infecteds 1
//                                         pop_mat(ind, 5) = I2t; // number of infecteds 2
//                                         pop_mat(ind, 6) = Rt;  // number of recovereds
//                                         
//                                         // increment row index
//                                         ind += 1;
//                                         
//                                         // check if any more transitions are possible (i.e. all rates equal to zero).
//                                         // if not, set keep_going to false.
//                                         if((rates[0] == 0) & (rates[1] == 0) & (rates[2] == 0)) {
//                                                   keep_going = false;
//                                         }
//                               }
//                     }
//           }
//           
//           return pop_mat;
// }
// 
// // Helper function for rbinding rows to pop_mat
// Rcpp::NumericMatrix rbind_popmat(Rcpp::NumericMatrix& M, int nrows) {
//           arma::mat M2 = Rcpp::as<arma::mat>(M);
//           M2.insert_rows(M2.n_rows, nrows);
//           return Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(M2));
// }
// 
// //' simulateSIR function, included in this file because of difficulty including sample.h in two separate files.
// //' Rcpp code to simulate a general stochastic epidemic and binomial counts.
// //'
// //' @param popsize population size
// //' @param obstimes vector of observation times
// //' @param params vector of parameters: beta, mu
// //' @param init_config initial configuration of compartment counts
// //'
// //' @return updated bookeeping objects
// // [[Rcpp::export]]
// Rcpp::NumericMatrix simulateSIRS(Rcpp::NumericVector obstimes, Rcpp::NumericVector params, Rcpp::IntegerVector init_config) {
//           
//           // initialize population bookkeeping matrix
//           int popsize = sum(init_config);
//           int n_events = 6 * popsize;
//           Rcpp::NumericMatrix pop_mat(n_events, 6);
//           pop_mat(0, 3) = init_config[0];
//           pop_mat(0, 4) = init_config[1];
//           pop_mat(0, 5) = init_config[2];
//           
//           // get parameters
//           double beta  = params[0];
//           double mu    = params[1];
//           double gamma = params[2];
//           
//           // set initial values
//           int nobs = obstimes.size();
//           double t = obstimes[0]; // current time
//           Rcpp::NumericVector dt(1); // time increment
//           double tmax = obstimes[nobs - 1]; // end of observation period
//           
//           double St = init_config[0]; // number of susceptibles
//           double It = init_config[1]; // number of infecteds
//           double Rt = init_config[2]; // number of recovereds
//           
//           // Vector of lumped rates
//           Rcpp::NumericVector probs(3);
//           Rcpp::NumericVector rates(3);
//           rates[0] = beta*St*It;
//           rates[1] = mu*It;
//           rates[2] = gamma*Rt;
//           
//           Rcpp::IntegerVector events = Rcpp::seq_len(2); // vector of event codes
//           Rcpp::IntegerVector next_event(1); // next event
//           Rcpp::IntegerVector subj_ind(1); // subject index
//           int subject(1); // subject ID
//           
//           Rcpp::IntegerVector susc_IDs = Rcpp::seq_len(St);            // IDs for susceptibles
//           Rcpp::IntegerVector infec_IDs = St + Rcpp::seq_len(It);      // IDs for infecteds
//           Rcpp::IntegerVector recov_IDs = St + It + Rcpp::seq_len(Rt); // IDs for recovereds
//           Rcpp::IntegerVector susc_Inds  = Rcpp::seq_len(St) - 1; // inds for susceptibles
//           Rcpp::IntegerVector infec_Inds = Rcpp::seq_len(It) - 1; // inds for infecteds
//           Rcpp::IntegerVector recov_Inds = Rcpp::seq_len(Rt) - 1; // inds for recovereds
//           std::vector<int> susceptibles = Rcpp::as<std::vector<int> >(susc_IDs); // convert to std vector
//           std::vector<int> infecteds = Rcpp::as<std::vector<int> >(infec_IDs);   // convert to std vector
//           std::vector<int> recovered = Rcpp::as<std::vector<int> >(recov_IDs);   // convert to std vector
//           std::vector<int> susc_inds = Rcpp::as<std::vector<int> >(susc_Inds);   // convert to std vector
//           std::vector<int> infec_inds = Rcpp::as<std::vector<int> >(infec_Inds); // convert to std vector
//           std::vector<int> recov_inds = Rcpp::as<std::vector<int> >(recov_Inds); // convert to std vector
//           
//           // set keep_going and the row index
//           bool keep_going = true;
//           int ind = 0;
//           
//           // start simulating
//           while(keep_going) {
//                     
//                     // sample the next time
//                     dt = Rcpp::rexp(1, sum(rates));
//                     t += dt[0];
//                     
//                     if(t > tmax) {
//                               keep_going = false; // stop simulating
//                               
//                     } else {
//                               // check if rows need to be added to the path matrix
//                               if(ind == n_events) {
//                                         pop_mat = rbind_popmat(pop_mat, n_events);
//                                         n_events = pop_mat.nrow();
//                               }
//                               
//                               // sample the next event
//                               probs = rates / sum(rates);
//                               next_event = Rcpp::RcppArmadillo::sample(events, 1, false, probs);
//                               
//                               if(next_event[0] == 1) { // infection occurs
//                                         
//                                         // choose subject
//                                         subj_ind = Rcpp::RcppArmadillo::sample(susc_inds, 1, false, Rcpp::NumericVector::create());
//                                         subject = susceptibles[subj_ind[0]];
//                                         
//                                         // move subject
//                                         susceptibles.erase(susceptibles.begin() + subj_ind[0]); // remove from susc IDs
//                                         infecteds.push_back(subject);                           // move to infecteds
//                                         
//                                         susc_inds.pop_back(); // shrink infected indices
//                                         infec_inds.push_back(infec_inds.size()); // expand infected indices
//                                         
//                                         // update compartment counts
//                                         St -= 1;
//                                         It += 1;
//                                         
//                                         // update rates
//                                         rates[0] = beta * St * It;
//                                         rates[1] = mu * It;
//                                         
//                                         // update bookkeeping matrix
//                                         pop_mat(ind, 0) = t; // new time
//                                         pop_mat(ind, 1) = subject; // subject ID
//                                         pop_mat(ind, 2) = 1; // event code
//                                         pop_mat(ind, 3) = St; // number of susceptibles
//                                         pop_mat(ind, 4) = It; // number of infecteds
//                                         pop_mat(ind, 5) = Rt; // number of recovereds
//                                         
//                                         // increment row index
//                                         ind += 1;
//                                         
//                                         // check if any more transitions are possible (i.e. all rates equal to zero).
//                                         // if not, set keep_going to false.
//                                         if((rates[0] == 0) & (rates[1] == 0) & (rates[2] == 0)) {
//                                                   keep_going = false;
//                                         }
//                                         
//                               } else if(next_event[0] == 2) { // recovery occurs
//                                         
//                                         // choose subject
//                                         subj_ind = Rcpp::RcppArmadillo::sample(infec_inds, 1, false, Rcpp::NumericVector::create());
//                                         subject = infecteds[subj_ind[0]];
//                                         
//                                         // move subject
//                                         infecteds.erase(infecteds.begin() + subj_ind[0]); // remove from infected IDs
//                                         recovered.push_back(subject);                     // move to recovered
//                                         
//                                         infec_inds.pop_back(); // shrink infected indices
//                                         recov_inds.push_back(recov_inds.size()); // expand infected indices
//                                         
//                                         // update compartment counts
//                                         It -= 1;
//                                         Rt += 1;
//                                         
//                                         // update rates
//                                         rates[0] = beta * St * It;
//                                         rates[1] = mu * It;
//                                         rates[2] = gamma * Rt;
//                                         
//                                         // update bookkeeping matrix
//                                         pop_mat(ind, 0) = t; // new time
//                                         pop_mat(ind, 1) = subject; // subject ID
//                                         pop_mat(ind, 2) = 2; // event code
//                                         pop_mat(ind, 3) = St; // number of susceptibles
//                                         pop_mat(ind, 4) = It; // number of infecteds
//                                         pop_mat(ind, 5) = Rt; // number of recovereds
//                                         
//                                         // increment row index
//                                         ind += 1;
//                                         
//                                         // check if any more transitions are possible (i.e. all rates equal to zero).
//                                         // if not, set keep_going to false.
//                                         if((rates[0] == 0) & (rates[1] == 0) & (rates[2] == 0)) {
//                                                   keep_going = false;
//                                         }
//                                         
//                               } else if(next_event[0] == 3) { // recovery occurs
//                                         
//                                         // choose subject
//                                         subj_ind = Rcpp::RcppArmadillo::sample(recov_inds, 1, false, Rcpp::NumericVector::create());
//                                         subject = recovered[subj_ind[0]];
//                                         
//                                         // move subject
//                                         recovered.erase(recovered.begin() + subj_ind[0]); // remove from recovered IDs
//                                         susceptibles.push_back(subject);                  // move to susceptibles
//                                         
//                                         recov_inds.pop_back(); // shrink infected indices
//                                         susc_inds.push_back(susc_inds.size()); // expand infected indices
//                                         
//                                         // update compartment counts
//                                         Rt -= 1;
//                                         St += 1;
//                                         
//                                         // update rates
//                                         rates[0] = beta * St * It;
//                                         rates[2] = gamma * Rt;
//                                         
//                                         // update bookkeeping matrix
//                                         pop_mat(ind, 0) = t; // new time
//                                         pop_mat(ind, 1) = subject; // subject ID
//                                         pop_mat(ind, 2) = 2; // event code
//                                         pop_mat(ind, 3) = St; // number of susceptibles
//                                         pop_mat(ind, 4) = It; // number of infecteds
//                                         pop_mat(ind, 5) = Rt; // number of recovereds
//                                         
//                                         // increment row index
//                                         ind += 1;
//                                         
//                                         // check if any more transitions are possible (i.e. all rates equal to zero).
//                                         // if not, set keep_going to false.
//                                         if((rates[0] == 0) & (rates[1] == 0) & (rates[2] == 0)) {
//                                                   keep_going = false;
//                                         }
//                               }
//                     }
//           }
//           
//           return pop_mat;
// }