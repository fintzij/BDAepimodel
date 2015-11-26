# this file contains pieces of code that are used in developing the package but
# are not necessary for the package to work.


# Initialize epimodel -----------------------------------------------------

set.seed(52787)
require(BDAepimodel)

popsize <- 750

r_meas_process <- function(state, meas_vars, params){
          rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
}

d_meas_process <- function(state, meas_vars, params, log = TRUE) {
          dbinom(x = state[, paste0(meas_vars, "_observed")], size = state[, paste0(meas_vars, "_augmented")], prob = params["rho"], log = log) 
}


epimodel <- init_epimodel(obstimes = seq(0, 10, by = 1),
                          popsize = popsize,
                          states = c("S", "I", "R"), 
                          params = c(beta = rnorm(1, 0.005, 1e-5), mu = rnorm(1, 1, 1e-5), rho = 0.5, S0 = 0.95, I0 = 0.01, R0 = 0.04), 
                          rates = c("beta * I", "mu"), 
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                          meas_vars = "I",
                          incidence_vars = "I",
                          r_meas_process = r_meas_process,
                          d_meas_process = d_meas_process)


# Simulate the epidemic ---------------------------------------------------

epimodel <- simulate_epimodel(epimodel = epimodel, lump = TRUE, trim = TRUE)

# Set the MCMC settings ---------------------------------------------------

# Kernel - MH - Beta, Mu, Rho ---------------------------------------------

# gibbs kernel
Rcpp::cppFunction("Rcpp::NumericVector getSuffStats(const Rcpp::NumericMatrix& pop_mat, const int ind_final_config) {
                  
                  // initialize sufficient statistics
                  int num_inf = 0;
                  int num_rec = 0;
                  double beta_suff = 0;
                  double mu_suff = 0;
                  
                  // initialize times
                  double cur_time = 0;
                  double next_time = 0;
                  
                  // compute the sufficient statistics
                  for(int j = 0; j < ind_final_config - 1; ++j) {
                  
                  cur_time = next_time;
                  next_time = pop_mat(j+1, 0); 
                  
                  beta_suff += pop_mat(j, 3) * pop_mat(j, 4) * (next_time - cur_time);
                  mu_suff += pop_mat(j, 4) * (next_time - cur_time);
                  
                  if(pop_mat(j + 1, 2) == 1) {
                  num_inf += 1;
                  }
                  
                  if(pop_mat(j + 1, 2) == 2) {
                  num_rec += 1;
                  }
                  }
                  
                  return Rcpp::NumericVector::create(num_inf, beta_suff, num_rec, mu_suff);
                  }")

gibbs_kernel <- function(epimodel) {
          
          # get sufficient statistics
          suff_stats          <- getSuffStats(epimodel$pop_mat, epimodel$ind_final_config)
          
          # update parameters
          proposal          <- epimodel$params
          proposal["beta"]  <- rgamma(1, 1 + suff_stats[1], 0.000001 + suff_stats[2])
          proposal["mu"]    <- rgamma(1, 1 + suff_stats[3], 0.000001 + suff_stats[4])
          proposal["rho"]   <- rbeta(1, shape1 = 1 + sum(epimodel$obs_mat[,"I_observed"]), shape2 = 1 + sum(epimodel$obs_mat[,"I_augmented"] - epimodel$obs_mat[,"I_observed"]))
          
          # update array of rate matrices
          epimodel            <- build_new_irms(epimodel, proposal)
          
          # update the eigen decompositions
          buildEigenArray(eigenvals = epimodel$eigen_values, eigenvecs = epimodel$eigen_vectors, inversevecs = epimodel$inv_eigen_vectors, irm_array = epimodel$irm)
          
          # get likelihoods under the new parameters
          pop_likelihood_new  <- calc_pop_likelihood(epimodel, log = TRUE)
          obs_likelihood_new  <- calc_obs_likelihood(epimodel, params = proposal, log = TRUE)
          
          # update parameters, likelihood objects, and eigen decompositions
          epimodel            <- update_params(epimodel, params = proposal, pop_likelihood = pop_likelihood_new, obs_likelihood = obs_likelihood_new)
          
          return(epimodel)
}


# save new parameters every iteration
# save every tenth configuration matrix
# resample 10 subject-level trajectories in between parameter updates
epimodel <- init_settings(epimodel,
                          niter = 100,
                          save_params_every = 1, 
                          save_configs_every = 5,
                          kernel = list(beta_mu_kernel, rho_kernel),
                          cov_mtx = diag(c(0.02, 0.02)),
                          configs_to_redraw = 70,
                          to_estimation_scale = to_estimation_scale,
                          from_estimation_scale = from_estimation_scale)

# epimodel <- init_tuning(epimodel,
#                         tune_params = c("beta", "mu", "rho"),
#                         cov_mtx = matrix(c(1, -0.865, -0.747, 
#                                            -0.865, 1, 0.831,
#                                            -0.747, 0.831, 1), nrow = 3),
#                         cov_scale = 0.02, 
#                         mintune = 20,
#                         maxtune = 50, 
#                         ntu = 5000,
#                         tunewt = 0.25,
#                         target_accept = c(0.15, 0.5))

# fit the model -----------------------------------------------------------

.epimodel <- epimodel
Rprof(prof <- tempfile())
.epimodel <- fit_epimodel(.epimodel, monitor = TRUE) 
Rprof()
summaryRprof(prof)

require(proftools)
plotProfileCallGraph(readProfileData(prof), score = "self")
