## ----load_packages, echo = F, include=F----------------------------------
library(BDAepimodel)
library(coda)
library(Rcpp)

## ----SIR_sim, warning=F, cache = F---------------------------------------
library(BDAepimodel)
library(Rcpp)

set.seed(81786)

obstimes  <- seq(1, 78, by=7)
params    <- c(beta = 0.0004, mu = 1/7, rho = 0.3, S0 = 0.948, I0 = 0.002, R0 = 0.05)
init_config <- c(1222,3,25)

epidemic <- simulateSIR(obstimes, params, init_config)
epidemic <- epidemic[epidemic[,1]!=0,]
epidemic <- rbind(c(1,0,0,init_config), epidemic, c(obstimes[obstimes>max(epidemic[,1])][1], 0, 0,epidemic[nrow(epidemic),4:6]))
obstimes <- obstimes[obstimes < max(epidemic[,1])]

dat <- cbind(obstimes, rbinom(length(obstimes), epidemic[findInterval(obstimes, epidemic[,1]),5], params["rho"]))
colnames(dat) <- c("time", "I")

plot(epidemic[,1], epidemic[,5], "l")
points(dat)


## ----initvals, echo = FALSE, results = 'asis'----------------------------
library(knitr)
inits <- t(data.frame(beta    = c(3.6, 3.8, 4, 4, 4.1, 4.2, 4.2, 6, 8.5) / 
                            c(1400, 1300, 1250, 1200, 1100, 900, 500, 300, 150) / 
                            c(7.25, 7.25, 7.5, 7.5, 7.75, 7.75, 8, 18, 25),
                    mu      = 1 / c(7.25, 7.25, 7.5, 7.5, 7.75, 7.75, 8, 15, 50)))
colnames(inits) = paste("N = ",c(1400, 1300, 1250, 1200, 1100, 900, 500, 300, 150), sep = "")
kable(inits, caption = "Initial rate parameter values.")

## ----SIR_kernel, warning = F, cache=F------------------------------------
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


gibbs_kernel <- function(epimodel) {
          
          # get sufficient statistics
          suff_stats          <- getSuffStats(epimodel$pop_mat, epimodel$ind_final_config)
          
          ### update parameters
          # beta ~ Gamma(0.00042 * 1250 / N, 1)
          # mu   ~ Gamma(0.35, 2)
          # rho  ~ Beta(1, 1)
          proposal          <- epimodel$params
          proposal["beta"]  <- rgamma(1, 0.00042 * (1250 / epimodel$popsize) + suff_stats[1], 1 + suff_stats[2])
          proposal["mu"]    <- rgamma(1, 0.35 + suff_stats[3], 2 + suff_stats[4])
          proposal["rho"]   <- rbeta(1, shape1 = 1 + sum(epimodel$obs_mat[,"I_observed"]), shape2 = 1 + sum(epimodel$obs_mat[,"I_augmented"] - epimodel$obs_mat[,"I_observed"]))
          
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
chain <- 1; popsize = 900 # both of these were varied by an external script
set.seed(52787 + chain)

# initialize the measurement process functions
r_meas_process <- function(state, meas_vars, params){
          
          rbinom(n = nrow(state), size = state[,meas_vars], prob = params["rho"])
          
}

d_meas_process <- function(state, meas_vars, params, log = TRUE) {
          dbinom(x = state[, "I_observed"], 
                 size = state[, "I_augmented"],
                 prob = params["rho"], log = log) 
}

# for setting the initial values
inits <- data.frame(popsize = c(1400, 1300, 1250, 1200, 1100, 900, 500, 300, 150),
                    beta    = c(3.6, 3.8, 4, 4, 4.1, 4.2, 4.2, 6, 8.5) / 
                            c(1400, 1300, 1250, 1200, 1100, 900, 500, 300, 150) / 
                            c(7.25, 7.25, 7.5, 7.5, 7.75, 7.75, 8, 18, 25),
                    mu      = 1 / c(7.25, 7.25, 7.5, 7.5, 7.75, 7.75, 8, 15, 50))

# initialize the epimodel object
epimodel <- init_epimodel(popsize = popsize,
                          states = c("S", "I", "R"), 
                          params = c(
                                    beta = rnorm(1, inits[inits$popsize == popsize, 2], 1e-5),
                                    mu = rnorm(1, inits[inits$popsize == popsize, 3], 1e-4),
                                    rho = rbeta(1,10,10),
                                    S0 = 0.97,
                                    I0 = 0.003,
                                    R0 = 0.027
                          ),
                          rates = c("beta * I", "mu"), 
                          flow = matrix(c(-1, 1, 0, 0, -1, 1), ncol = 3, byrow = T), 
                          dat = dat,
                          time_var = "time",
                          meas_vars = "I",
                          initdist_prior = c(100, 1, 5),
                          r_meas_process = r_meas_process,
                          d_meas_process = d_meas_process)

# initialize the model object
epimodel <- init_settings(epimodel,
                          niter = 10, # set to 100,000 in the paper
                          save_params_every = 1, 
                          save_configs_every = 5000,
                          kernel = list(gibbs_kernel),
                          configs_to_redraw = ceiling(0.1 * popsize),
                          analytic_eigen = "SIR",
                          ecctmc_method = "unif",
                          52787 + chain)

# fit the model
epimodel <- fit_epimodel(epimodel, monitor = TRUE)


