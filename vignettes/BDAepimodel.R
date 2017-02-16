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
                 size = state[,meas_vars], # binomial sample of the unobserved prevalenc
                 prob = params["rho"])     # sampling probability
}

d_meas_process <- function(state, meas_vars, params, log = TRUE) {
          # note that the names of the measurement variables are endowed with suffixes "_observed" and "_augmented". This is required.
          # we will declare the names of the measurement variables shortly.
          dbinom(x = state[, "I_observed"], 
                 size = state[, "I_augmented"], 
                 prob = params["rho"], log = log)
}

## ------------------------------------------------------------------------
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
Rcpp::cppFunction("Rcpp::NumericVector getSuffStats_SIR(const Rcpp::NumericMatrix& pop_mat, const int ind_final_config) {
                  
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
gibbs_kernel_SIR <- function(epimodel) {
          
          # get sufficient statistics using the previously compiled getSuffStats_SIR function (above)
          suff_stats <- getSuffStats_SIR(epimodel$pop_mat, epimodel$ind_final_config)
          
          # update parameters from their univariate full conditional distributions
          # Priors: beta ~ gamma(0.3, 1000)
          #         mu   ~ gamma(1, 8)
          #         rho  ~ beta(2, 7)
          proposal          <- epimodel$params # params is the vector of ALL model parameters
          proposal["beta"]  <- rgamma(1, 0.3 + suff_stats[1], 1000 + suff_stats[2])
          proposal["mu"]    <- rgamma(1, 1 + suff_stats[3], 8 + suff_stats[4])
          proposal["rho"]   <- rbeta(1, shape1 = 2 + sum(epimodel$obs_mat[, "I_observed"]),
                                        shape2 = 7 + sum(epimodel$obs_mat[, "I_augmented"] - epimodel$obs_mat[, "I_observed"]))
          
          # update array of rate matrices
          epimodel <- build_new_irms(epimodel, proposal)
          
          # update the eigen decompositions (This function is built in and computes eigen decompositions analytically)
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

## ------------------------------------------------------------------------
# grab the data that was simulated previously. No need to redefine the measurement process functions, they remain unchanged.
dat <- epimodel$dat 

chain <- 1 # this was set by a batch script that ran chains 1, 2, and 3 in parallel
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
                          meas_vars = "I",                                                     # name of measurement variable
                          initdist_prior = c(90, 2, 5), ### Parameters for the dirichlet prior distribution for the initial state probs
                          r_meas_process = r_meas_process,
                          d_meas_process = d_meas_process)

## ---- warning = FALSE----------------------------------------------------
epimodel <- init_settings(epimodel,
                          niter = 10,  # this was set to 100,000 in the paper
                          save_params_every = 1, 
                          save_configs_every = 2, # this was set to 250 in the paper 
                          kernel = list(gibbs_kernel_SIR),
                          configs_to_redraw = 1, # this was set to 75 in the paper
                          analytic_eigen = "SIR", # compute eigen decompositions and matrix inverses analytically
                          ecctmc_method = "unif")   # sample subject paths in interevent intervals via modified rejection sampling

epimodel <- fit_epimodel(epimodel, monitor = FALSE)


## ------------------------------------------------------------------------
head(epimodel$results$params)
ts.plot(epimodel$results$params[,"beta"], ylab = expression(beta))
plot(hist(epimodel$results$params[,"rho"]), main = "Binomial sampling probability")

## ---- fig.cap=""---------------------------------------------------------
plot_latent_posterior(epimodel, 
                      states = "I", times = epimodel$obstimes, cm = "mn")

