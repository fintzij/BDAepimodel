## ----echo=T,eval=T-------------------------------------------------------
library(BDAepimodel)
library(coda)
library(Rcpp)
library(pomp)
library(gtools)
library(ggplot2)

## ----echo=T,eval=T-------------------------------------------------------
data<-c(7,2,4,14,14,20,21,8,7,6,1,2,2,0,0,1)

## ----echo=T,eval=T-------------------------------------------------------
rmeas<-"
cases=rbinom(I,exp(rho)/(1+exp(rho)));//represent the data
"
dmeas<-"
lik=dbinom(cases,I,exp(rho)/(1+exp(rho)),give_log); //return the loglikelihood
"

## ----echo=T,eval=T-------------------------------------------------------
sir.step<-"
double rate[2];
double dN[2];
rate[0]=exp(beta)*I;    //Infection rate
rate[1]=exp(mu);       //recovery rate
reulermultinom(1,S,&rate[0],dt,&dN[0]);   //generate the number of newly infected people
reulermultinom(1,I,&rate[1],dt,&dN[1]);   //generate the number of newly recovered people
if(!R_FINITE(S)) Rprintf(\"%lg %lg %lg %lg %lg %lg %lg %lg %lg\\n\",dN[0],rate[0],dN[1],rate[1],beta,mu,S,I,R);
S+=-dN[0];          //update the number of Susceptible
I+=dN[0]-dN[1];     //update the number of Infection
R+=dN[1];           //update the number of Recovery
" 

## ----echo=T,eval=T-------------------------------------------------------
sir <- pomp(
  data = data.frame(cases = data, time = seq(1, 106, by = 7)),  #"cases" is the dataset, "time" is the observation time
  times = "time",
  t0 = 1,  #initial time point
  dmeasure = Csnippet(dmeas),  
  rmeasure = Csnippet(rmeas),  #return the likelihood
  rprocess = euler.sim(step.fun = Csnippet(sir.step), delta.t = 1/12), #delta.t is here
  statenames = c("S", "I", "R"),  #state space variable name
  paramnames = c("beta", "mu", "rho", "theta1", "theta2"), #parameters name
  initializer = function(params, t0, ...) {
    ps <-exp(params["theta1"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  # initial prob of susceptible
    pi <- 1 / (1 + exp(params["theta1"]) + exp(params["theta2"]))                     # initial prob of infection
    pr <-exp(params["theta2"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  # initial prob of recovery
    return(setNames(as.numeric(rmultinom(
      1, 750, prob = c(ps, pi, 1 - ps - pi)
    )), c("S", "I", "R")))
  },
  params = c(    #for data generation(We didn't generate data here, but pomp need this part.)
    beta = -5,  
    mu = -2,
    rho = 0,
    theta1 = 6,
    theta2 = 1
  )
)

## ----echo=T,eval=T-------------------------------------------------------
sir.dprior <- function(params, ..., log) {
          f <- dgamma(exp(params[1]), shape = 0.3, rate = 1000, log = TRUE) + params[1] + # log prior for log(infection rate) "beta"
                    dgamma(exp(params[2]), shape = 1, rate = 8, log = TRUE) + params[2] +    # log prior for log(recovery rate) "mu"
                    dbeta(exp(params[3]) / (1 + exp(params[3])), 2, 7, log = TRUE) + 
                    params[3] - log((1 + exp(params[3])) ^ 2) + #log prior for logit(sampling probablity) "rho"
                    log(ddirichlet(c(exp(params[4]) / (1 + exp(params[4]) + exp(params[5])),
                                     1 / (1 + exp(params[4]) + exp(params[5])),
                                     exp(params[5]) / (1 + exp(params[4]) + exp(params[5]))), c(90, 2, 5))) + 
                    params[4] + params[5] - 3 * log(1 + exp(params[4]) + exp(params[5])) #log prior for logit(initial value)
          
          if (log) {
                    f         
          } else{
                    exp(f)
          }
} 


## ----echo=T,eval=T-------------------------------------------------------
param.initial <- c(
          beta = log(rnorm(1, 0.00035, 1e-5)),
          mu = log(rnorm(1, 1/7, 1e-2)),
          rho = -1,
          theta1 = 4,
          theta2 = 1
)

## ----echo=T,eval=T,warning=F---------------------------------------------
pmcmc1 <- pmcmc(
          pomp(sir, dprior = sir.dprior),
          #given the prior function
          start = param.initial,
          #given the initial value of the parameters
          Nmcmc = 10,
          #number of mcmc steps
          Np = 100,
          max.fail = Inf,
          proposal = mvn.rw.adaptive(0.1 * c(
                    beta = 0.2,
                    #sampling variance for beta
                    mu = 0.2,
                    #sampling variance for mu
                    rho = 0.3,
                    #sampling variance for rho
                    theta1 = 0.2,
                    #sampling variance for theta1
                    theta2 = 0.5 #sampling variance for theta2
          ), 
          scale.start = 100, 
          shape.start = 100)
)

## ----echo=T,eval=T,warning=F---------------------------------------------
start_time <- Sys.time();  #calculation of time
pmcmc1 <- pmcmc(
          pmcmc1,
          #given the prior function
          start = param.initial,
          #given the initial value of the parameters
          Nmcmc = 10,
          max.fail = Inf,
          proposal = mvn.rw(covmat(pmcmc1))
)

end_time <- Sys.time(); #calculation of time
run_time <- difftime(end_time, start_time, units = "hours") #calculation of time

pomp_results <- list(time = run_time, results = pmcmc1) 

## ----echo=T,eval=T-------------------------------------------------------
rmeas<-"
cases=rbinom(I,exp(rho)/(1+exp(rho)));//represent the data
"
dmeas<-"
lik=dbinom(cases,I,exp(rho)/(1+exp(rho)),give_log); //return the loglikelihood
"

## ----echo=T,eval=T-------------------------------------------------------
# define the rate function, stoichiometry matrix, and rate-event dependency matrix
SIR_stoich <- cbind(infection = c(-1, 1, 0), # infections yield S-1, I+1
                    recovery  = c(0, -1, 1)) # recoveries yield I-1, R+1

SIR_depmat <- cbind(infection = c(1, 1, 0), # infectivity rate updated by changes to S and I
                    recovery  = c(0, 1, 0)) # recovery rate updated by changes to I
# SIR rate function
# j the number of the elementary event (1 = infection, 2 = recovery)
# x named numeric vector with the value of the state of the process at time t
# t time
# params named numeric vector containing the parameters
# returns single numerical value with the rate of the elementary event
SIR_rates <- function(j, x, t, params, ...) {
          switch(j,
                 exp(params["beta"]) * x["S"] * x["I"], # infection
                 exp(params["mu"]) * x["I"] # recovery
          )
}
# instatiate the gillespie stepper function
SIR_sim <- gillespie.sim(rate.fun = SIR_rates,
                         v = SIR_stoich,
                         d = SIR_depmat)

## ----echo=T,eval=T-------------------------------------------------------
sir <- pomp(
  data = data.frame(time = seq(1, 106, by = 7), cases = data),  #"cases" is the dataset, "time" is the observation time
  times = "time",
  t0 = 1,                      # initial time point
  dmeasure = Csnippet(dmeas),  # evaluates the density of the measurement process
  rmeasure = Csnippet(rmeas),  # simulates from the measurement process
  rprocess = SIR_sim,          # simulates from the latent process
  statenames = c("S", "I", "R"),  #state space variable name
  paramnames = c("beta", "mu", "rho", "theta1", "theta2"), #parameters name
  initializer = function(params, t0, ...) {
            ps <-exp(params["theta1"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  # initial prob of susceptible
            pi <- 1 / (1 + exp(params["theta1"]) + exp(params["theta2"]))                     # initial prob of infection
            pr <-exp(params["theta2"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  # initial prob of recovery
            return(setNames(as.numeric(rmultinom(1, 750, prob = c(ps, pi, 1 - ps - pi))), c("S", "I", "R")))
  },
  params = c(    #for data generation(We didn't generate data here, but pomp need this part.)
    beta = -5,  
    mu   = -2,
    rho  = 0,
    theta1 = 6,
    theta2 = 1
  )
)

## ----echo=T,eval=T-------------------------------------------------------
sir.dprior <- function(params, ..., log) {
  f <- dgamma(exp(params[1]), shape = 0.3, rate = 1000, log = TRUE) + params[1] + # log prior for log(infection rate) "beta"
            dgamma(exp(params[2]), shape = 1, rate = 8, log = TRUE) + params[2] +    # log prior for log(recovery rate) "mu"
            dbeta(exp(params[3]) / (1 + exp(params[3])), 2, 7, log = TRUE) + 
            params[3] - log((1 + exp(params[3])) ^ 2) + #log prior for logit(sampling probablity) "rho"
            log(ddirichlet(c(exp(params[4]) / (1 + exp(params[4]) + exp(params[5])),
                             1 / (1 + exp(params[4]) + exp(params[5])),
                             exp(params[5]) / (1 + exp(params[4]) + exp(params[5]))), c(90, 2, 5))) + 
            params[4] + params[5] - 3 * log(1 + exp(params[4]) + exp(params[5])) #log prior for logit(initial value)
    
          if (log) {
                    f         
          } else{
                    exp(f)
          }
} 

## ----echo=T,eval=T-------------------------------------------------------
param.initial <- c(
  beta = log(rnorm(1, 0.00035, 1e-5)),
  mu = log(rnorm(1, 1/7, 1e-2)),
  rho = -1,
  theta1 = 4,
  theta2 = 1
)

## ----echo=T,eval=T,warning=F---------------------------------------------
pmcmc1 <- pmcmc(
          pomp(sir, dprior = sir.dprior),
          #given the prior function
          start = param.initial,
          #given the initial value of the parameters
          Nmcmc = 10,
          #number of mcmc steps
          Np = 200,
          max.fail = Inf,
          proposal = mvn.rw.adaptive(0.1 * c(
                    beta = 0.2,
                    #sampling variance for beta
                    mu = 0.2,
                    #sampling variance for mu
                    rho = 0.3,
                    #sampling variance for rho
                    theta1 = 0.2,
                    #sampling variance for theta1
                    theta2 = 0.5 #sampling variance for theta2
          ), 
          scale.start = 100, 
          shape.start = 100)
)

## ----echo=T,eval=T,warning=F---------------------------------------------
start_time <- Sys.time();  #calculation of time
pmcmc1 <- pmcmc(
          pmcmc1,
          #given the prior function
          start = param.initial,
          #given the initial value of the parameters
          Nmcmc = 10,
          max.fail = Inf,
          proposal = mvn.rw(covmat(pmcmc1))
)

end_time <- Sys.time(); #calculation of time
run_time <- difftime(end_time, start_time, units = "hours") #calculation of time

pomp_results <- list(time = run_time, results = pmcmc1) 

