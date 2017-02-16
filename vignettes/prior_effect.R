## ----load_packages, echo = F, include=F----------------------------------
library(BDAepimodel)
library(coda)
library(Rcpp)

## ----SIR_sim, warning=F, cache = F---------------------------------------
set.seed(52787)

# declare the functions for simulating from and evaluating the log-density of the measurement process
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
          dbinom(x = state[, "I_observed"], 
                 size = state[, "I_augmented"], 
                 prob = params["rho"], log = log)
}

# initialize the stochastic epidemic model object
epimodel <- init_epimodel(obstimes = seq(0, 105, by = 7),                             # vector of observation times
                          popsize = 750,                                              # population size
                          states = c("S", "I", "R"),                                  # compartment names
                          params = c(beta = 0.00035,                                  # infectivity parameter
                                     mu = 1/7,                                        # recovery rate
                                     rho = 0.2,                                       # binomial sampling probability
                                     S0 = 0.9, I0 = 0.03, R0 = 0.07),                 # initial state probabilities
                          rates = c("beta * I", "mu"),                                # unlumped transition rates
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),  # flow matrix
                          meas_vars = "I",                                            # name of measurement variable
                          r_meas_process = r_meas_process,                            # measurement process functions
                          d_meas_process = d_meas_process)

# simulate the epidemic and the dataset.  
epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)

dat <- epimodel$dat 
true_path <- epimodel$pop_mat

plot(x = epimodel$pop_mat[,"time"], y = epimodel$pop_mat[,"I"], xlim = c(0,85), "l", xlab = "Time", ylab = "Prevalence")
points(x = epimodel$dat[,"time"], y = epimodel$dat[,"I"])

## ----SIR_kernel, warning = F, cache=F------------------------------------
# define the hyperprior parameters for the rates and sampling probability
beta_prior <- matrix(c(3, 10000, 0.3, 1000), nrow = 2); colnames(beta_prior) <- c("informative", "diffuse")
mu_prior   <- matrix(c(3, 20, 0.1, 0.8), nrow = 2); colnames(mu_prior) <- c("informative", "diffuse")
rho_prior  <- matrix(c(21, 75, 1, 1), nrow = 2); colnames(rho_prior) <- c("informative", "diffuse") 


# set the prior for this vignette - these were set with an external batch script
rates_prior <- 1; # "informative"
samp_prior  <- 1; # "informative"

# helper function for computing the sufficient statistics for the SIR model rate parameters
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

# MCMC transition kernel for the SIR model rate parameters and the binomial
# sampling probability. The prior distributions for the parameters are contained
# in this function.

gibbs_kernel <- function(epimodel) {
          
          # get sufficient statistics using the previously compiled getSuffStats function (above)
          suff_stats          <- getSuffStats(epimodel$pop_mat, epimodel$ind_final_config)
          
          # update parameters from their univariate full conditional distributions
          # Priors: beta ~ gamma(0.3, 1000)
          #         mu   ~ gamma(1, 8)
          #         rho  ~ beta(21, 75)
          proposal          <- epimodel$params # params is the vector of ALL model parameters
          proposal["beta"]  <- rgamma(1, beta_prior[1,rates_prior] + suff_stats[1], beta_prior[2,rates_prior] + suff_stats[2])
          proposal["mu"]    <- rgamma(1, mu_prior[1,rates_prior] + suff_stats[3], mu_prior[2,rates_prior] + suff_stats[4])
          proposal["rho"]   <- rbeta(1, 
                                     shape1 = rho_prior[1,samp_prior] + sum(epimodel$obs_mat[, "I_observed"]),
                                     shape2 = rho_prior[2,samp_prior] + sum(epimodel$obs_mat[, "I_augmented"] - 
                                                                                    epimodel$obs_mat[, "I_observed"]))
          
          # update array of rate matrices
          epimodel          <- build_new_irms(epimodel, proposal)
          
          # update the eigen decompositions (This function is built in)
          buildEigenArray_SIR(real_eigenvals = epimodel$real_eigen_values,
                              imag_eigenvals = epimodel$imag_eigen_values,
                              eigenvecs      = epimodel$eigen_vectors, 
                              inversevecs    = epimodel$inv_eigen_vectors, 
                              irm_array      = epimodel$irm, 
                              n_real_eigs    = epimodel$n_real_eigs, 
                              initial_calc   = FALSE)
          
          # get log-likelihood of the observations under the new parameters
          obs_likelihood_new  <- calc_obs_likelihood(epimodel, params = proposal, log = TRUE) #### NOTE - log = TRUE
          
          # get the new population level CTMC log-likelihood
          pop_likelihood_new  <- epimodel$likelihoods$pop_likelihood_cur +
                    suff_stats[1] * (log(proposal["beta"]) - log(epimodel$params["beta"])) +
                    suff_stats[3] * (log(proposal["mu"]) - log(epimodel$params["mu"])) -
                    suff_stats[2] * (proposal["beta"] - epimodel$params["beta"]) - 
                    suff_stats[4] * (proposal["mu"] - epimodel$params["mu"])
          
          # update parameters, likelihood objects, and eigen decompositions
          epimodel  <-
                    update_params(
                              epimodel,
                              params = proposal,
                              pop_likelihood = pop_likelihood_new,
                              obs_likelihood = obs_likelihood_new
                    )
          
          return(epimodel)
}

## ----SIR_inference, warning=F, cache=F, messages = F---------------------
chain <- 1 # set by an external script

set.seed(52787 + chain)

# initial values for initial state parameters
init_dist <- MCMCpack::rdirichlet(1, c(9,0.5,0.1))
epimodel <- init_epimodel(popsize = 750,                                                       # population size
                          states = c("S", "I", "R"),                                           # compartment names
                          params = c(beta = abs(rnorm(1, 0.00035, 5e-5)),                      # per-contact infectivity rate
                                     mu = abs(rnorm(1, 1/7, 0.02)),                            # recovery rate
                                     rho = rbeta(1, 21, 75),                                   # binomial sampling probability
                                     S0 = init_dist[1], I0 = init_dist[2], R0 = init_dist[3]), # initial state probabilities
                          rates = c("beta * I", "mu"),                                         # unlumped transition rates
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),           # flow matrix
                          dat = dat,                                                           # dataset
                          time_var = "time",                                                   # name of time variable in the dataset
                          meas_vars = "I",                                                     # name of measurement var in the dataset
                          initdist_prior = c(9,0.2,0.5), ### Parameters for the dirichlet prior distribution for the initial state probs
                          r_meas_process = r_meas_process,
                          d_meas_process = d_meas_process)

epimodel <- init_settings(epimodel,
                          niter = 10, # set to 100000 for the paper
                          save_params_every = 1, 
                          save_configs_every = 250, 
                          kernel = list(gibbs_kernel),
                          configs_to_redraw = 5, # set to 75 for the paper
                          analytic_eigen = "SIR")

epimodel <- fit_epimodel(epimodel, monitor = FALSE)


