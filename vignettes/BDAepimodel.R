## ---- eval = FALSE-------------------------------------------------------
#  library(devtools)
#  install_github("fintzij/BDAepimodel",build_vignettes=TRUE)
#  library(BDAepimodel)
#  

## ---- include=FALSE------------------------------------------------------
require(ggplot2)
require(Rcpp)
require(MCMCpack)
require(BDAepimodel)
set.seed(1834)

## ------------------------------------------------------------------------
r_meas_process <- function(state, meas_vars, params){
          # in our example, rho will be the name of the binomial sampling probability parameter.
          # this function returns a matrix of observed counts
          rbinom(n = nrow(state), 
                 size = state[,meas_vars],
                 prob = params["rho"])
}

d_meas_process <- function(state, meas_vars, params, log = TRUE) {
          # note that the names of the measurement variables are endowed with suffixes "_observed" and "_augmented". This is required.
          # we will declare the names of the measurement variables shortly.
          dbinom(x = state[, paste0(meas_vars, "_observed")], 
                 size = state[, paste0(meas_vars, "_augmented")], 
                 prob = params["rho"], log = log) 
}

## ------------------------------------------------------------------------
epimodel <- init_epimodel(obstimes = seq(0, 15, by = 1),                              # vector of observation times
                          popsize = 750,                                              # population size
                          states = c("S", "I", "R"),                                  # compartment names
                          params = c(beta = 2.5/(750 * 0.9) * 0.5,                    # infectivity parameter
                                     mu = 0.5,                                        # recovery rate
                                     rho = 0.2,                                       # binomial sampling probability
                                     S0 = 0.9, I0 = 0.03, R0 = 0.07),                 # initial state probabilities
                          rates = c("beta * I", "mu"),                                # transition rates
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),  # flow matrix
                          meas_vars = "I",                                            # name of measurement variable
                          r_meas_process = r_meas_process,                            # measurement process functions
                          d_meas_process = d_meas_process)

## ------------------------------------------------------------------------
# simulate the epidemic and the dataset. The lump argument instructs the simulation function to simulate the epidemic on the lumped state space of counts, which is more efficient than doing so on the unlumped state space of individuals. 
epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)
plot(x = epimodel$pop_mat[,"time"], y = epimodel$pop_mat[,"I"], "l", ylim = c(-5, 200), xlab = "Time", ylab = "Prevalence")
points(x = epimodel$dat[,"time"], y = epimodel$dat[,"I"])

## ------------------------------------------------------------------------
head(epimodel$obs_mat)

## ------------------------------------------------------------------------
head(epimodel$init_config)

## ------------------------------------------------------------------------
head(epimodel$pop_mat)

## ------------------------------------------------------------------------
Rcpp::cppFunction("Rcpp::NumericVector getSuffStats(const Rcpp::NumericMatrix& pop_mat, const int ind_final_config) {
                  
          // initialize sufficient statistics
          int num_inf = 0;       // number of infection events
          int num_rec = 0;       // number of recovery events
          double beta_suff = 0;  // integrated hazard for the infectivity
          double mu_suff = 0;    // integrated hazard for the recovery

          // initialize times
          double cur_time = 0;              // current time
          double next_time = pop_mat(0,0);  // time of the first event
          double dt = 0;                    // time increment
          
          // compute the sufficient statistics - loop through the pop_mat matrix until
          // reaching the row for the final observation time
          for(int j = 0; j < ind_final_config - 1; ++j) {
          
                    cur_time = next_time;         
                    next_time = pop_mat(j+1, 0); // grab the time of the next event
                    dt = next_time - cur_time;   // compute the time increment
                    
                    beta_suff += pop_mat(j, 3) * pop_mat(j, 4) * dt; // add S*I*(t_{j+1} - t_j) to beta_suff
                    mu_suff += pop_mat(j, 4) * dt;                   // add I*(t_{j+1} - t_j) to mu_suff
                    
                    // increment the count for the next event
                    if(pop_mat(j + 1, 2) == 1) {  
                              num_inf += 1;
                    } else if(pop_mat(j + 1, 2) == 2) {
                              num_rec += 1;
                    }
          }
                  
          // return the vector of sufficient statistics for the rate parameters
          return Rcpp::NumericVector::create(num_inf, beta_suff, num_rec, mu_suff);
}")

## ------------------------------------------------------------------------
gibbs_kernel <- function(epimodel) {
          
          # get sufficient statistics using the previously compiled getSuffStats function (above)
          suff_stats          <- getSuffStats(epimodel$pop_mat, epimodel$ind_final_config)
          
          # update parameters from their univariate full conditional distributions
          # Priors: beta ~ gamma(0.002, 1)
          #         mu   ~ gamma(2.1, 4)
          #         rho  ~ beta(2/3, 1)
          proposal          <- epimodel$params # params is the vector of ALL model parameters
          proposal["beta"]  <- rgamma(1, 0.002 + suff_stats[1], 1 + suff_stats[2])
          proposal["mu"]    <- rgamma(1, 2.1 + suff_stats[3], 4 + suff_stats[4])
          proposal["rho"]   <- rbeta(1, shape1 = 2/3 + sum(epimodel$obs_mat[,"I_observed"]), shape2 = 1 + sum(epimodel$obs_mat[,"I_augmented"] - epimodel$obs_mat[,"I_observed"]))
          
          
          #### THE FOLLOWING SECTION CAN BE INCLUDED PIECEMEAL AS APPROPRIATE
          # update array of rate matrices - MUST BE INCLUDED WHEN UPDATING RATE PARAMETERS
          epimodel            <- build_new_irms(epimodel, proposal)
          
          # update the eigen decompositions - MUST BE INCLUDED WHEN UPDATING RATE PARAMETERS
          buildEigenArray_SIR(real_eigenvals = epimodel$real_eigen_values,
                              imag_eigenvals = epimodel$imag_eigen_values,
                              eigenvecs      = epimodel$eigen_vectors, 
                              inversevecs    = epimodel$inv_eigen_vectors, 
                              irm_array      = epimodel$irm, 
                              n_real_eigs    = epimodel$n_real_eigs, 
                              initial_calc   = FALSE)          
          # get log-likelihoods under the new parameters
          # MUST BE INCLUDED WHEN UPDATING RATE PARAMETERS
          pop_likelihood_new  <- calc_pop_likelihood(epimodel, log = TRUE) #### NOTE - log = TRUE
          
          # MUST BE INCLUDED WHEN UPDATING SAMPLING PARAMETERS
          obs_likelihood_new  <- calc_obs_likelihood(epimodel, params = proposal, log = TRUE) #### NOTE - log = TRUE
          
          # update parameters, likelihood objects, and eigen decompositions
          epimodel <-
                    update_params(
                    epimodel,
                    params = proposal,
                    pop_likelihood = pop_likelihood_new,
                    obs_likelihood = obs_likelihood_new
                    )
                    
          return(epimodel)
}


## ------------------------------------------------------------------------
# grab the data that was simulated previously. No need to redefine the measurement process functions, they remain unchanged.
dat <- epimodel$dat 

# initial values for the initial state probability parameters
init_dist <- MCMCpack::rdirichlet(1, c(9,0.5,0.5))

# initialize the epimodel object
epimodel <- init_epimodel(popsize = 750,                                                       # population size
                          states = c("S", "I", "R"),                                           # compartment names
                          params = c(beta = rnorm(1, 2.8 / 750, 1e-3),                         # per-contact infectivity rate
                                     mu = rnorm(1, 0.5, 0.1),                                  # recovery rate
                                     rho = 150 / 750,                                          # binomial sampling probability
                                     S0 = init_dist[1], I0 = init_dist[2], R0 = init_dist[3]), # initial state probabilities
                          rates = c("beta * I", "mu"),                                         # unlumped transition rates
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),           # flow matrix
                          dat = dat,                                                           # dataset
                          time_var = "time",                                                   # name of time variable in the dataset
                          meas_vars = "I",                                                     # name of measurement var in the dataset
                          initdist_prior = c(9,0.2,0.5), ### Parameters for the dirichlet prior distribution for the initial state probs
                          r_meas_process = r_meas_process,
                          d_meas_process = d_meas_process)


