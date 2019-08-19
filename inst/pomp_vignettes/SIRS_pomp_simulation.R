## ----echo=T,eval=T-------------------------------------------------------
library(BDAepimodel)
library(coda)
library(Rcpp)
library(pomp)
library(gtools)
library(ggplot2)

## ----echo=T,eval=T-------------------------------------------------------
data<-c(1,4,6,6,4,9,12,23,29,32,29,48,32,22,15,13,8,7,3,2,1,0,0,1,3,1,0,1,3,3,3,1,5,3,5,3,2,1,0,6,2,3,6,3,6,5,5,8,15,17,14,8,8)

## ----echo=T,eval=T-------------------------------------------------------
rmeas<-"
cases=rbinom(I,exp(rho)/(1+exp(rho)));//represent the data
"
dmeas<-"
lik=dbinom(cases,I,exp(rho)/(1+exp(rho)),give_log); //return the loglikelihood
"

## ----echo=T,eval=T-------------------------------------------------------
sirs.step<-"
double rate[3];
double dN[3];
rate[0]=exp(beta)*I;  //Infection rate
rate[1]=exp(mu);     //Recovery rate
rate[2]=exp(gamma);    //Loss of immunity rate
reulermultinom(1,S,&rate[0],dt,&dN[0]); //generate the number of newly infected people
reulermultinom(1,I,&rate[1],dt,&dN[1]); //generate the number of newly infected people
reulermultinom(1,R,&rate[2],dt,&dN[2]); //generate the number of people who newly lose immunity
if(!R_FINITE(S)) Rprintf(\"%lg %lg %lg %lg %lg %lg %lg %lg %lg\\n\",dN[0],rate[0],dN[1],rate[1],beta,mu,S,I,R);
S+=-dN[0]+dN[2];  //update the number of Susceptible
I+=dN[0]-dN[1];   //update the number of Infection
R+=dN[1]-dN[2];   //update the number of Recovery
"

## ----echo=T,eval=T-------------------------------------------------------
sirs <- pomp(
  data = data.frame(cases = data, time = seq(1, 365, by = 7)), #"cases" is the dataset, "time" is the observation time
  times = "time",
  t0 = 1,  #initial time point
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),  #return the likelihood
  rprocess = euler.sim(step.fun = Csnippet(sirs.step), delta.t = 1/3), #delta.t is here
  statenames = c("S", "I", "R"),  #state space variable name
  paramnames = c("beta", "mu", "gamma", "rho", "theta1", "theta2"), #parameters name
  initializer = function(params, t0, ...) {
    ps <-exp(params["theta1"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  #given the initial proportion of susceptible
    pi <- 1 / (1 + exp(params["theta1"]) + exp(params["theta2"]))                     #given the initial proportion of infection
    pr <-exp(params["theta2"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  #given the initial proportion of recovery
    return(setNames(as.numeric(rmultinom(1, 200, prob = c(ps, pi, pr))), c("S", "I", "R")))
  },
  # not used, but pomp needs this
  params = c(
            beta = log(0.0001),
            mu = log(1 / 14),
            gamma = log(1 / 140),
            rho = 1.5,
            theta1 = log(49),
            theta2 = log(0.00000001)
  )
)

## ----echo=T,eval=T-------------------------------------------------------
sirs.dprior <- function(params, ..., log) {
          f <- (
                    dgamma(exp(params[1]), 0.1, 100, log = TRUE) + params[1] + #log prior for log(infection rate) "beta"
                    dgamma(exp(params[2]), 1.8, 14, log = TRUE) + params[2] + #log prior for log(recovery rate) "mu"
                    dgamma(exp(params[3]), 0.0625, 10, log = TRUE) + params[3] + #log prior for log(loss of susceptible) "gamma"
                    dbeta(exp(params[4]) / (1 + exp(params[4])), 5, 1, log = TRUE) + 
                    params[4] - log((1 + exp(params[4])) ^ 2) + #log prior for logit(sampling probablity) "rho"
                    log(ddirichlet(c(exp(params[5]) / (1 + exp(params[5]) + exp(params[6])),
                                     1 / (1 + exp(params[5]) + exp(params[6])),
                                     exp(params[6]) / (1 + exp(params[5]) + exp(params[6]))),
                                   c(9, 0.15, 0.001))) +
                    params[5] + params[6] - 3 * log(1 + exp(params[5]) + exp(params[6])) #log prior for logit(initial value)
          )
          if (log) {
                    f
          } else {
                    exp(f)
          }
}

## ----echo=T,eval=T-------------------------------------------------------
param.initial <- c(
          beta = log(abs(rnorm(1, 0.00077, 1e-4))),
          mu = log(abs(rnorm(1, 1/18, 1e-4))),
          gamma = log(abs(rnorm(1, 1/140, 1e-4))),
          rho = 1.5,
          theta1 = 3,
          theta2 = -5
)

## ----echo=T,eval=T,warning=F---------------------------------------------
pmcmc1 <- pmcmc(
          pomp(sirs, dprior = sirs.dprior),
          #given the prior function
          start = param.initial,
          #given the initial value of the parameters
          Nmcmc = 10,
          #number of mcmc steps
          Np = 500,
          max.fail = Inf,
          proposal = mvn.rw.adaptive(0.1 * c(
                    beta = 0.1,  #sampling variance for beta
                    mu = 0.1,    #sampling variance for mu
                    gamma = 0.1, #sampling variance for gamma
                    rho = 0.1,   #sampling variance for rho
                    theta1 = 0.1,#sampling variance for theta1
                    theta2 = 0.1 #sampling variance for theta2
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
SIRS_stoich <- cbind(infection = c(-1, 1, 0), # infections yield S-1, I+1
                     recovery  = c(0, -1, 1), # recoveries yield I-1, R+1
                     susceptibility = c(1, 0, -1)) # susceptibility events yield R-1, S+1

SIRS_depmat <- cbind(infection = c(1, 1, 0), # infectivity rate updated by changes to S and I
                     recovery  = c(0, 1, 0), # recovery rate updated by changes to I
                     susceptibility = c(0, 0, 1)) # susceptibility rate updated by changes to R 
# SIRS rate function
# j the number of the elementary event (1 = infection, 2 = recovery, 3 = loss of immunity)
# x named numeric vector with the value of the state of the process at time t
# t time
# params named numeric vector containing the parameters
# returns single numerical value with the rate of the elementary event
SIRS_rates <- function(j, x, t, params, ...) {
          switch(j,
                 exp(params["beta"]) * x["S"] * x["I"], # infection
                 exp(params["mu"]) * x["I"], # recovery
                 exp(params["gamma"]) * x["R"] # loss of immunity
          )
}

# instatiate the gillespie stepper function
SIRS_sim <- gillespie.sim(rate.fun = SIRS_rates,
                          v = SIRS_stoich,
                          d = SIRS_depmat)

## ----echo=T,eval=T-------------------------------------------------------
# initialize the pomp object
sirs <- pomp(
          data = data.frame(time = seq(1, 365, by = 7), cases = data),  #"cases" is the dataset, "time" is the observation time
          times = "time",
          t0 = 1,                      # initial time point
          dmeasure = Csnippet(dmeas),  # evaluates the density of the measurement process
          rmeasure = Csnippet(rmeas),  # simulates from the measurement process
          rprocess = SIRS_sim,          # simulates from the latent process
          statenames = c("S", "I", "R"),  #state space variable name
          paramnames = c("beta", "mu", "gamma", "rho", "theta1", "theta2"), #parameters name
          initializer = function(params, t0, ...) {
                    ps <-exp(params["theta1"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  # initial prob of susceptible
                    pi <- 1 / (1 + exp(params["theta1"]) + exp(params["theta2"]))                     # initial prob of infection
                    pr <-exp(params["theta2"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  # initial prob of recovery
                    return(setNames(as.numeric(rmultinom(1, 200, prob = c(ps, pi, pr))), c("S", "I", "R")))
          },
          # not used, but pomp needs this
          params = c(
                    beta = log(0.0001),
                    mu = log(1 / 14),
                    gamma = log(1 / 150),
                    rho = log(0.85),
                    theta1 = log(49),
                    theta2 = log(0.00000001)
          )
)

## ----echo=T,eval=T-------------------------------------------------------
sirs.dprior <- function(params, ..., log) {
          f <- (
          dgamma(exp(params[1]), 0.1, 100, log = TRUE) + params[1] + #log prior for log(infection rate) "beta"
          dgamma(exp(params[2]), 1.8, 14, log = TRUE) + params[2] + #log prior for log(recovery rate) "mu"
          dgamma(exp(params[3]), 0.0625, 10, log = TRUE) + params[3] + #log prior for log(loss of susceptible) "gamma"
          dbeta(exp(params[4]) / (1 + exp(params[4])), 5, 1, log = TRUE) + 
              params[4] - log((1 + exp(params[4])) ^ 2) + #log prior for logit(sampling probablity) "rho"
          log(ddirichlet(c(exp(params[5]) / (1 + exp(params[5]) + exp(params[6])),
                       1 / (1 + exp(params[5]) + exp(params[6])),
                       exp(params[6]) / (1 + exp(params[5]) + exp(params[6]))),
                     c(90, 1.5, 0.01))) +
          params[5] + params[6] - 3 * log(1 + exp(params[5]) + exp(params[6])) #log prior for logit(initial value)
          )
          if (log) {
                    f
          } else {
                    exp(f)
          }
}

## ----echo=T,eval=T-------------------------------------------------------
param.initial <- c(
  beta = log(abs(rnorm(1, 0.00077, 1e-4))),
  mu = log(abs(rnorm(1, 1/18, 1e-4))),
  gamma = log(abs(rnorm(1, 1/140, 1e-4))),
  rho = 1.5,
  theta1 = 3,
  theta2 = -5
)

## ----echo=T,eval=T,warning=F---------------------------------------------
pmcmc1 <- pmcmc(
          pomp(sirs, dprior = sirs.dprior),
          #given the prior function
          start = param.initial,
          #given the initial value of the parameters
          Nmcmc = 10,
          #number of mcmc steps
          Np = 500,
          max.fail = Inf,
          proposal = mvn.rw.adaptive(0.1 * c(
                    beta = 0.1,
                    #sampling variance for beta
                    mu = 0.1,
                    #sampling variance for mu
                    gamma = 0.1,
                    #sampling variance for gamma
                    rho = 0.1,
                    #sampling variance for rho
                    theta1 = 0.1,
                    #sampling variance for theta1
                    theta2 = 0.1   #sampling variance for theta2
          ), 
          scale.start = 100, 
          shape.start = 100, 
          max.scaling = 1.1)
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

