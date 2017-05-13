## ----load_packages, echo = F, include=F----------------------------------
library(BDAepimodel)
library(coda)
library(MASS)
library(ggplot2)
library(Rcpp)

## ----data, echo = FALSE, results = 'asis'--------------------------------
set.seed(52787)
st.date <- as.Date("1978/1/22")
en.date <- as.Date("1978/2/4")
date.seq <- seq(st.date,en.date, by = "day")

bbs <- data.frame(time = rep(date.seq,1),
                  count = c(rep("Observed count", 14)),
                  I = c(3,8,28,76,222,293,257,237,192,126,70,28,12,5))

theme_set(theme_bw())
ggplot(bbs, aes(x = time, y = I, colour = count, shape = count)) + geom_point(size = 4) + geom_line(size = 1) + labs(x = "Date (1978)", y = "Number of infected boys") + theme(text = element_text(size = 20)) + scale_colour_discrete(guide = "none") + scale_shape(guide = "none")

# The fitting function requires a matrix, whereas ggplot requires a data frame.
# We're just creating a matrix version of the data, no big deal.
dat <- matrix(c(0:13, 3,8,28,76,222,293,257,237,192,126,70,28,12,5), ncol = 2)
colnames(dat) <- c("time", "I")

## ----meas_procs----------------------------------------------------------
# binomial emissions
d_meas_process_binom <- function(state, meas_vars, params, log = TRUE) {
          dbinom(x = state[, "I_observed"], 
                 size = state[, "I_augmented"],
                 prob = params["rho"], log = log)
}

# negative binomial emissions
d_meas_process_negbinom <- function(state, meas_vars, params, log = TRUE) {
          dnbinom(x = state[, "I_observed"], 
                  size = params["phi"], 
                  mu = state[, "I_augmented"] * params["rho"], log = log)
}

r_meas_process_binom <- function(state, meas_vars, params){
          rbinom(n = nrow(state),
                 size = state[, meas_vars],
                 prob = params["rho"])
}

r_meas_process_negbinom <- function(state, meas_vars, params){
          rnbinom(n = nrow(state), 
                  size = params["phi"],
                  mu = state[,meas_vars] * params["rho"])
}


## ----init_epimodels------------------------------------------------------
SIR_binom <- init_epimodel(popsize = 763,
                          states = c("S", "I", "R"),
                          params = c(
                                    beta = rnorm(1, 0.002, 1e-4),
                                    mu = rnorm(1, 0.35, 1e-3),
                                    rho = rbeta(1, 100, 10),
                                    S0 = 0.99,
                                    I0 = 0.005,
                                    R0 = 0.005),
                          rates = c("beta * I", "mu"),
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),
                          dat = dat,
                          time_var = "time",
                          meas_vars = "I",
                          initdist_prior = c(900, 3, 9),
                          r_meas_process = r_meas_process_binom,
                          d_meas_process = d_meas_process_binom)

SEIR_binom <- init_epimodel(popsize = 763,
                          states = c("S", "E", "I", "R"),
                          params = c(
                                    beta = rnorm(1, 0.002, 1e-4),
                                    gamma = rnorm(1, 5, 1e-1),
                                    mu = rnorm(1, 0.35, 1e-3),
                                    rho = rbeta(1, 100, 10),
                                    S0 = 0.99,
                                    E0 = 0.007,
                                    I0 = 0.005,
                                    R0 = 0.005),
                          rates = c("beta * I", "gamma", "mu"),
                          flow = matrix(c(-1, 1, 0, 0, 
                                          0, -1, 1, 0, 
                                          0, 0, -1, 1), ncol = 4, byrow = T),
                          dat = dat,
                          time_var = "time",
                          meas_vars = "I",
                          initdist_prior = c(900, 6, 3, 9),
                          r_meas_process = r_meas_process_binom,
                          d_meas_process = d_meas_process_binom)

SIR_negbinom <- init_epimodel(popsize = 763,
                          states = c("S", "I", "R"),
                          params = c(
                                    beta = rnorm(1, 0.002, 1e-4),
                                    mu = rnorm(1, 0.35, 1e-3),
                                    rho = rbeta(1, 25, 0.3),
                                    phi = runif(1, 1,10),
                                    S0 = 0.99,
                                    I0 = 0.005,
                                    R0 = 0.005),
                          rates = c("beta * I", "mu"),
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T),
                          dat = dat,
                          time_var = "time",
                          meas_vars = "I",
                          initdist_prior = c(900, 3, 9),
                          r_meas_process = r_meas_process_negbinom,
                          d_meas_process = d_meas_process_negbinom)


SEIR_negbinom <- init_epimodel(popsize = 763,
                          states = c("S", "E", "I", "R"),
                          params = c(
                                    beta = rnorm(1, 0.002, 1e-4),
                                    gamma = rnorm(1, 5, 1e-1),
                                    mu = rnorm(1, 0.35, 1e-3),
                                    rho = rbeta(1, 100, 10),
                                    phi = runif(1, 1,10),
                                    S0 = 0.99,
                                    E0 = 0.007,
                                    I0 = 0.005,
                                    R0 = 0.005),
                          rates = c("beta * I", "gamma", "mu"),
                          flow = matrix(c(-1, 1, 0, 0, 
                                          0, -1, 1, 0, 
                                          0, 0, -1, 1), ncol = 4, byrow = T),
                          dat = dat,
                          time_var = "time",
                          meas_vars = "I",
                          initdist_prior = c(900, 6, 3, 9),
                          r_meas_process = r_meas_process_negbinom,
                          d_meas_process = d_meas_process_negbinom)


## ----SIR_kernel, warning = F, cache=F------------------------------------
# helper function for computing the sufficient statistics for the SIR model rate parameters
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

### Helper function for computing the SEIR model sufficient statistics
Rcpp::cppFunction("Rcpp::NumericVector getSuffStats_SEIR(const Rcpp::NumericMatrix& pop_mat, const int ind_final_config) {
                  
          // initialize sufficient statistics
          int num_exp = 0;       // number of exposure events
          int num_inf = 0;       // number of exposed --> infectious events
          int num_rec = 0;       // number of recovery events
          double beta_suff  = 0; // integrated hazard for the exposure
          double gamma_suff = 0; // integrated hazard for addition of infectives
          double mu_suff    = 0; // integrated hazard for the recovery
          
          // initialize times
          double cur_time = 0;              // current time
          double next_time = pop_mat(0,0);  // time of the first event
          double dt = 0;                    // time increment

          // compute the sufficient statistics - loop through the pop_mat matrix until
          // reaching the row for the final observation time
          for(int j = 0; j < ind_final_config - 1; ++j) {

                    cur_time = next_time;         
                    next_time = pop_mat(j+1, 0); // grab the time of the next event
                    dt = next_time - cur_time;

                    beta_suff  += pop_mat(j, 3) * pop_mat(j, 5) * dt; // add S*I*(t_{j+1} - t_j) to beta_suff
                    gamma_suff += pop_mat(j, 4) * dt;                 // add E*(t_{j+1} - t_j) to gamma_suff
                    mu_suff    += pop_mat(j, 5) * dt;                 // add I*(t_{j+1} - t_j) to mu_suff
                    
                    if(pop_mat(j + 1, 2) == 1) {  
                              num_exp += 1;            // if the next event is an exposure, increment the number of exposures
                    } else if(pop_mat(j + 1, 2) == 2) {
                              num_inf += 1;            // if the next event adds an infective, increment the number of infections
                    } else if(pop_mat(j + 1, 2) == 3) {
                              num_rec += 1;            // if the next event is a recover, increment the number of recovery
                    }
          }
          
          // return the vector of sufficient statistics for the rate parameters
          return Rcpp::NumericVector::create(num_exp, beta_suff, num_inf, gamma_suff, num_rec, mu_suff);
}")

# MCMC transition kernel for the SIR model rate parameters and the binomial
# sampling probability. The prior distributions for the parameters are contained
# in this function.
SIR_binom_kernel <- function(epimodel) {

          # get sufficient statistics
          suff_stats          <- getSuffStats_SIR(epimodel$pop_mat, epimodel$ind_final_config)

          # update parameters
          proposal          <- epimodel$params
          proposal["beta"]  <- rgamma(1, 0.001 + suff_stats[1], 1 + suff_stats[2])
          proposal["mu"]    <- rgamma(1, 1 + suff_stats[3], 2 + suff_stats[4])
          proposal["rho"]   <- rbeta(1, shape1 = 2 + sum(epimodel$obs_mat[, "I_observed"]),
                                        shape2 = 1 + sum(epimodel$obs_mat[, "I_augmented"] - epimodel$obs_mat[, "I_observed"]))
          
          # update array of rate matrices
          epimodel            <- build_new_irms(epimodel, proposal)

          # update the eigen decompositions
          buildEigenArray_SIR(real_eigenvals = epimodel$real_eigen_values,
                              imag_eigenvals = epimodel$imag_eigen_values,
                              eigenvecs      = epimodel$eigen_vectors, 
                              inversevecs    = epimodel$inv_eigen_vectors, 
                              irm_array      = epimodel$irm, 
                              n_real_eigs    = epimodel$n_real_eigs, 
                              initial_calc   = FALSE)
          
          # get likelihoods under the new parameters
          # get log-likelihood of the observations under the new parameters
          obs_likelihood_new  <- calc_obs_likelihood(epimodel, params = proposal, log = TRUE) #### NOTE - log = TRUE
          
          # get the new population level CTMC log-likelihood
          pop_likelihood_new  <- epimodel$likelihoods$pop_likelihood_cur +
                    suff_stats[1] * (log(proposal["beta"]) - log(epimodel$params["beta"])) +
                    suff_stats[3] * (log(proposal["mu"]) - log(epimodel$params["mu"])) -
                    suff_stats[2] * (proposal["beta"] - epimodel$params["beta"]) - 
                    suff_stats[4] * (proposal["mu"] - epimodel$params["mu"])
          
          # update parameters, likelihood objects, and eigen decompositions
          epimodel <- update_params(
                              epimodel,
                              params = proposal,
                              pop_likelihood = pop_likelihood_new,
                              obs_likelihood = obs_likelihood_new
                    )
          
          return(epimodel)
}

# MCMC transition kernel for the SEIR model rate parameters and the binomial
# sampling probability. The prior distributions for the parameters are contained
# in this function.
SEIR_binom_kernel <- function(epimodel) {
          
          # get sufficient statistics using the previously compiled getSuffStats function (above)
          suff_stats <- getSuffStats_SEIR(epimodel$pop_mat, epimodel$ind_final_config)
          
          # update parameters from their univariate full conditional distributions
          proposal          <- epimodel$params # params is the vector of ALL model parameters
          proposal["beta"]  <- rgamma(1, 0.001 + suff_stats[1], 1 + suff_stats[2])
          proposal["gamma"] <- rgamma(1, 0.001 + suff_stats[3], 1 + suff_stats[4])
          proposal["mu"]    <- rgamma(1, 1 + suff_stats[5], 2 + suff_stats[6])
          proposal["rho"]   <- rbeta(1, 
                                     shape1 = 2 + sum(epimodel$obs_mat[,"I_observed"]), 
                                     shape2 = 1 + sum(epimodel$obs_mat[,"I_augmented"]- epimodel$obs_mat[,"I_observed"]))
          
          # update array of rate matrices
          epimodel          <- build_new_irms(epimodel, proposal)
          
          # update the eigen decompositions
          buildEigenArray_SEIR(real_eigenvals = epimodel$real_eigen_values,
                               imag_eigenvals = epimodel$imag_eigen_values,
                               eigenvecs      = epimodel$eigen_vectors, 
                               inversevecs    = epimodel$inv_eigen_vectors, 
                               irm_array      = epimodel$irm, 
                               n_real_eigs    = epimodel$n_real_eigs, 
                               initial_calc   = FALSE)
          
          # get the data log-likelihood under the new parameters
          obs_likelihood_new <- calc_obs_likelihood(epimodel, params = proposal, log = TRUE) #### NOTE - log = TRUE
          
          # compute the new population level CTMC log-likelihood
          pop_likelihood_new <- epimodel$likelihoods$pop_likelihood_cur +
                    suff_stats[1] * (log(proposal["beta"]) - log(epimodel$params["beta"])) + 
                    suff_stats[3] * (log(proposal["gamma"]) - log(epimodel$params["gamma"])) +
                    suff_stats[5] * (log(proposal["mu"]) - log(epimodel$params["mu"])) -
                    suff_stats[2] * (proposal["beta"] - epimodel$params["beta"]) - 
                    suff_stats[4] * (proposal["gamma"] - epimodel$params["gamma"]) - 
                    suff_stats[6] * (proposal["mu"] - epimodel$params["mu"])
          
          # update parameters, likelihood objects, and eigen decompositions
          epimodel <- update_params(epimodel,
                                    params = proposal,
                                    pop_likelihood = pop_likelihood_new,
                                    obs_likelihood = obs_likelihood_new
          )
          
          return(epimodel)
}

# MCMC transition kernel for the SIR model rate parameters. The prior distributions for the parameters are contained
# in this function.
SIR_negbinom_rates_kernel <- function(epimodel) {
          
          # get sufficient statistics using the previously compiled getSuffStats function (above)
          suff_stats          <- getSuffStats_SIR(epimodel$pop_mat, epimodel$ind_final_config)
          
          # update parameters from their univariate full conditional distributions
          # Priors: beta ~ gamma(0.003, 10)
          #         mu   ~ gamma(1, 8)
          proposal          <- epimodel$params # params is the vector of ALL model parameters
          proposal["beta"]  <- rgamma(1, 0.001 + suff_stats[1], 1 + suff_stats[2])
          proposal["mu"]    <- rgamma(1, 1 + suff_stats[3], 2 + suff_stats[4])
          
          # update array of rate matrices
          epimodel <- build_new_irms(epimodel, proposal)
          
          # update the eigen decompositions
          buildEigenArray_SIR(real_eigenvals = epimodel$real_eigen_values,
                              imag_eigenvals = epimodel$imag_eigen_values,
                              eigenvecs      = epimodel$eigen_vectors, 
                              inversevecs    = epimodel$inv_eigen_vectors, 
                              irm_array      = epimodel$irm, 
                              n_real_eigs    = epimodel$n_real_eigs, 
                              initial_calc   = FALSE)
          
          # get the new population level CTMC log-likelihood
          pop_likelihood_new  <- epimodel$likelihoods$pop_likelihood_cur +
                    suff_stats[1] * (log(proposal["beta"]) - log(epimodel$params["beta"])) +
                    suff_stats[3] * (log(proposal["mu"]) - log(epimodel$params["mu"])) -
                    suff_stats[2] * (proposal["beta"] - epimodel$params["beta"]) - 
                    suff_stats[4] * (proposal["mu"] - epimodel$params["mu"])
          
          # update parameters, likelihood objects, and eigen decompositions
          epimodel <- update_params(epimodel, params = proposal, pop_likelihood = pop_likelihood_new)
          
          return(epimodel)
}

# MCMC transition kernel for the SIR model rate parameters. The prior distributions for the parameters are contained
# in this function.
SEIR_negbinom_rates_kernel <- function(epimodel) {
          
          # get sufficient statistics using the previously compiled getSuffStats function (above)
          suff_stats          <- getSuffStats_SEIR(epimodel$pop_mat, epimodel$ind_final_config)
          
          # update parameters from their univariate full conditional distributions
          # Priors: beta ~ gamma(0.003, 10)
          #         mu   ~ gamma(1, 8)
          proposal          <- epimodel$params # params is the vector of ALL model parameters
          proposal["beta"]  <- rgamma(1, 0.001 + suff_stats[1], 1 + suff_stats[2])
          proposal["gamma"] <- rgamma(1, 0.001 + suff_stats[3], 1 + suff_stats[4])
          proposal["mu"]    <- rgamma(1, 1 + suff_stats[5], 2 + suff_stats[6])
          
          # update array of rate matrices
          epimodel          <- build_new_irms(epimodel, proposal)
          
          # update the eigen decompositions
          buildEigenArray_SEIR(real_eigenvals = epimodel$real_eigen_values,
                               imag_eigenvals = epimodel$imag_eigen_values,
                               eigenvecs      = epimodel$eigen_vectors, 
                               inversevecs    = epimodel$inv_eigen_vectors, 
                               irm_array      = epimodel$irm, 
                               n_real_eigs    = epimodel$n_real_eigs, 
                               initial_calc   = FALSE)
          
          # get the new population level CTMC log-likelihood
          pop_likelihood_new <- epimodel$likelihoods$pop_likelihood_cur +
                    suff_stats[1] * (log(proposal["beta"]) - log(epimodel$params["beta"])) + 
                    suff_stats[3] * (log(proposal["gamma"]) - log(epimodel$params["gamma"])) +
                    suff_stats[5] * (log(proposal["mu"]) - log(epimodel$params["mu"])) -
                    suff_stats[2] * (proposal["beta"] - epimodel$params["beta"]) - 
                    suff_stats[4] * (proposal["gamma"] - epimodel$params["gamma"]) - 
                    suff_stats[6] * (proposal["mu"] - epimodel$params["mu"])
          
          # update parameters, likelihood objects, and eigen decompositions
          epimodel <- update_params(epimodel, params = proposal, pop_likelihood = pop_likelihood_new)
          
          return(epimodel)
}

## RWMH transition ernels for the negative binomial parameters
# functions for converting rho and phi to and from the estimation scale for the negative binomial parameters
to_estimation_scale <- function(params, epimodel) {
          params["rho"] <- logit(params["rho"])
          params["phi"] <- log(params["phi"])
          return(params)
}

from_estimation_scale <- function(params_est, epimodel) {
          params_est["rho"] <- expit(params_est["rho"])
          params_est["phi"] <- exp(params_est["phi"])
          return(params_est)
}

RWMH_kernel <- function(epimodel) {
          
          # get existing parameter values - proposal will contain the parameters on the estimation scale
          # while proposal_nat contains the proposal on the natural scale
          proposal_est <- proposal_nat <- epimodel$params
          
          # set log likelihood function
          log_likelihood_cur <- epimodel$likelihoods$obs_likelihood
          
          # make proposal - edit this line appropriately
          proposal_est[c("rho", "phi")] <- epimodel$sim_settings$to_estimation_scale(proposal_nat[c("rho", "phi")], epimodel)
          
          # compute the prior probability of the current parameter values
          prior_prob_cur <- sum(
                    dbeta(proposal_nat["rho"], 2, 1, log = TRUE) - proposal_est["rho"] - 2*log(1 + exp(-proposal_est["rho"])),
                    dgamma(proposal_nat["phi"], 1, 0.1, log = TRUE) + proposal_est["phi"]
          )
          
          # make proposal
          proposal_est[c("rho", "phi")] <- MASS::mvrnorm(n = 1, mu = proposal_est[c("rho", "phi")], Sigma = epimodel$sim_settings$cov_mtx)
          
          # transform back to the natural scale
          proposal_nat[c("rho", "phi")] <- epimodel$sim_settings$from_estimation_scale(proposal_est[c("rho", "phi")], epimodel)
          
          # get prior likelihood for proposed parameters - edit appropriately
          prior_prob_new <- sum(
                    dbeta(proposal_nat["rho"], 2, 1, log = TRUE) - proposal_est["rho"] - 2*log(1 + exp(-proposal_est["rho"])),
                    dgamma(proposal_nat["phi"], 1, 0.1, log = TRUE) + proposal_est["phi"]
          )
          
          # get population level complete data likelihood for the new parameters
          # if the population-level CTMC likelihood or the measurement process
          # is unchanged, comment out the corresponding line
          obs_likelihood_new  <- calc_obs_likelihood(epimodel = epimodel, params = proposal_nat, log = TRUE)
          log_likelihood_new  <- obs_likelihood_new
          
          # compute the acceptance probability and accept or reject proposal 
          accept_prob <- log_likelihood_new - log_likelihood_cur + prior_prob_new - prior_prob_cur
          
          if(accept_prob >= 0 || accept_prob >= log(runif(1))) {
                    
                    # update parameters and likelihood objects
                    epimodel <-
                              update_params(
                                        epimodel,
                                        params = proposal_nat,
                                        obs_likelihood = obs_likelihood_new
                              )
          } 
          
          return(epimodel)
}

## ----inference, warning=F, cache=F, messages = F-------------------------
chain <- 1 # this was set by a batch script that ran chains 1, 2, and 3 in parallel
set.seed(52787 + chain)


SIR_binom <- init_settings(SIR_binom,
                          niter = 10, # this was set to 100000 for the paper
                          save_params_every = 1,
                          save_configs_every = 250,
                          kernel = list(SIR_binom_kernel),
                          configs_to_redraw = 10, # this was set to 100 for the paper
                          ecctmc_method = "unif",
                          analytic_eigen = "SIR")

SEIR_binom <- init_settings(SEIR_binom,
                          niter = 10, # this was set to 100000 for the paper
                          save_params_every = 1,
                          save_configs_every = 250,
                          kernel = list(SEIR_binom_kernel),
                          configs_to_redraw = 10, # this was set to 100 for the paper
                          ecctmc_method = "unif",
                          analytic_eigen = "SEIR")

SIR_negbinom <- init_settings(SIR_negbinom,
                          niter = 10, # this was set to 100000 for the paper
                          save_params_every = 1,
                          save_configs_every = 250,
                          kernel = list(SIR_negbinom_rates_kernel, RWMH_kernel), 
                          cov_mtx = 0.5 * matrix(c(5.779382, -1.0413166,
                                             -1.0413166,  0.7541288), 2),
                          configs_to_redraw = 10, # this was set to 100 for the ppaer 
                          to_estimation_scale = to_estimation_scale,
                          from_estimation_scale = from_estimation_scale,
                          ecctmc_method = "unif",
                          analytic_eigen = "SIR")

SEIR_negbinom <- init_settings(SEIR_negbinom,
                          niter = 10, # this was set to 100000 for the paper
                          save_params_every = 1,
                          save_configs_every = 250,
                          kernel = list(SEIR_negbinom_rates_kernel, RWMH_kernel), 
                          cov_mtx = 0.5 * matrix(c(5.779382, -1.0413166,
                                             -1.0413166,  0.7541288), 2),
                          configs_to_redraw = 10, # this was set to 100000 for the paper 
                          to_estimation_scale = to_estimation_scale,
                          from_estimation_scale = from_estimation_scale,
                          ecctmc_method = "unif",
                          analytic_eigen = "SEIR")

# Fit the models --------------------------------------------------------

SIR_binom_results     <- fit_epimodel(SIR_binom, monitor = TRUE)
SEIR_binom_results    <- fit_epimodel(SEIR_binom, monitor = TRUE)
SIR_negbinom_results  <- fit_epimodel(SIR_negbinom, monitor = TRUE)
SEIR_negbinom_results <- fit_epimodel(SEIR_negbinom, monitor = TRUE)


