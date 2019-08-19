## ----echo=T,eval=T-------------------------------------------------------
library(BDAepimodel)
library(coda)
library(Rcpp)
library(pomp)
library(gtools)
library(ggplot2)

## ----echo=T,eval=T-------------------------------------------------------
set.seed(1834)
data<-c(0,0,0,1,0,1,1,1,2,1,2,1,0,3,1,0,2,0,2,1,2,2,2,1,0,2,1,0,0,2,0,1,1,0,1,1,1,0,1,2,2,1,1,2,4,2,4,3,5,3,1,1,5,5,4,3,1,1,3,2,0,0,2,2,2,3,2,5,4,1,3,1,4,3,1,2,2,5,2,4,2,1,2,1,3,1,1,3,1,1,0,2,2,3,5,1,0,1,0,1,0,0,0,0,1)

## ----echo=T,eval=T-------------------------------------------------------
rmeas<-"
cases=rbinom(I,exp(rho)/(1+exp(rho)));//represent the data
"
dmeas<-"
lik=dbinom(cases,I,exp(rho)/(1+exp(rho)),give_log); //return the loglikelihood
"

## ----echo=T,eval=T-------------------------------------------------------
seir.step<-"
double rate[3];
double dN[3];
rate[0]=exp(beta)*I;   //Infection
rate[1]=exp(gamma);    //Rate of S to E
rate[2]=exp(mu);       //Rate of E to I
reulermultinom(1,S,&rate[0],dt,&dN[0]);   //generate the number of newly from S to E
reulermultinom(1,E,&rate[1],dt,&dN[1]);   //generate the number of newly from E to I
reulermultinom(1,I,&rate[2],dt,&dN[2]);   //generate the number of newly from I to R
if(!R_FINITE(S)) Rprintf(\"%lg %lg %lg %lg %lg %lg %lg %lg %lg\\n\",dN[0],rate[0],dN[1],rate[1],beta,mu,S,I,R);
S+=-dN[0];       //update the number of Susceptible
E+=dN[0]-dN[1];  //update the number of E
I+=dN[1]-dN[2];  //update the number of I
R+=dN[2];        //update the number of R
"

## ----echo=T,eval=T,warning=F---------------------------------------------
seir <- pomp(
  data = data.frame(cases = data, time = seq(1, 729, by = 7)), #"cases" is the dataset, "time" is the observation time
  times = "time",
  t0 = 1,  #initial time point
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
  rprocess = euler.sim(step.fun = Csnippet(seir.step), delta.t = 1/12), # tau-leaping over 2 hour intervals
  statenames = c("S", "E", "I", "R"), #state space variable name
  paramnames = c("beta", "gamma", "mu", "rho", "theta1", "theta2", "theta3"), #parameters name
  initializer = function(params, t0, ...) {
            #Initial proportion of S, E, I, R
            ps <- exp(params["theta1"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
            pe <- exp(params["theta2"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
            pi <- 1 / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
            pr <- exp(params["theta3"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))  
          return(setNames(as.numeric(rmultinom(1, 500, prob = c(ps, pe, pi, pr))), c("S", "E", "I", "R")))
  },
  params = c(
    beta = log(0.00072),
    gamma = log(0.5),
    mu = log(0.2),
    rho = log(1 / 4),
    theta1 = log(90),
    theta2 = log(2),
    theta3 = log(7)
  )
)

## ----echo=T,eval=T-------------------------------------------------------
trans<-function(a,b,c){
  a_t<-exp(a)
  b_t<-exp(b)
  c_t<-exp(c)
  return(a_t*b_t*c_t/((1+a_t+b_t+c_t)^4))
}


seir.dprior <- function(params, ..., log) {
          f <- (dgamma(exp(params[1]), 1, 10000, log = TRUE) + params[1] +  #log prior for beta
                          dgamma(exp(params[2]), 1, 11, log = TRUE) + params[2] +  #log prior for gamma
                          dgamma(exp(params[3]), 3.2, 100, log = TRUE) + params[3] +   #log prior for mu
                          dbeta(exp(params[4]) / (1 + exp(params[4])), 3.5, 6.5, log = TRUE) +
                          params[4] - log((1 + exp(params[4])) ^ 2) +      #log prior for rho
                          log(ddirichlet(c(exp(params[5]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                                           exp(params[6]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                                           1 / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                                           exp(params[7]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7]))),
                                         c(100, 0.01, 0.4, 0.01)) * trans(params[5], params[6], params[7]))  #log prior for initial SEIR value
          )
          if (log) {
                    f
          } else {
                    exp(f)
          }
}

## ----echo=T,eval=T-------------------------------------------------------
param.initial<- c(
  beta = log(abs(rnorm(1, 0.00005, 1e-6))),
  gamma = log(abs(rnorm(1, 0.15, 0.01))),
  mu = log(abs(rnorm(1, 0.04, 0.001))),
  rho = -1,
  theta1 = 5,
  theta2 = 1,
  theta3 = 1
)

## ----echo=T,eval=T,warning=F---------------------------------------------
pmcmc1 <- pmcmc(
          pomp(seir, dprior = seir.dprior),
          #given the prior function
          start = param.initial,
          #given the initial value of the parameters
          Nmcmc = 10,
          #number of mcmc steps
          Np = 200,
          max.fail = Inf,
          proposal = mvn.rw.adaptive(
                    0.5 * c(
                              beta = 0.1,
                              #sampling variance for beta
                              gamma = 0.1,
                              #sampling variance for gamma
                              mu = 0.1,
                              #sampling variance for mu
                              rho = 0.1,
                              #sampling variance for rho
                              theta1 = 0.1,
                              #sampling variance for theta1
                              theta2 = 0.1,
                              #sampling variance for theta2
                              theta3 = 0.1 #sampling variance for theta3
                    ), 
                    scale.start = 100, 
                    shape.start = 100
          )
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
SEIR_stoich <- cbind(exposure  = c(-1, 1, 0, 0), # exposures yield S-1, E+1
                     infection = c(0, -1, 1, 0), # infections yield E-1, I+1
                     recovery  = c(0, 0, -1, 1)) # recoveries yield I-1, R+1

SEIR_depmat <- cbind(exposure  = c(1, 0, 1, 0), # exposure rate updated by changes to S and I
                     infection = c(0, 1, 0, 0), # infection rate updated by changes to E
                     recovery  = c(0, 0, 1, 0)) # recovery rate updated by changes to I
# SIR rate function
# j the number of the elementary event (1 = infection, 2 = recovery)
# x named numeric vector with the value of the state of the process at time t
# t time
# params named numeric vector containing the parameters
# returns single numerical value with the rate of the elementary event
SEIR_rates <- function(j, x, t, params, ...) {
          switch(j,
                 exp(params["beta"]) * x["S"] * x["I"], # exposure
                 exp(params["gamma"]) * x["E"], # infection
                 exp(params["mu"]) * x["I"] # recovery
          )
}

# instatiate the gillespie stepper function
SEIR_sim <- gillespie.sim(rate.fun = SEIR_rates,
                          v = SEIR_stoich,
                          d = SEIR_depmat)

## ----echo=T,eval=T,warning=F---------------------------------------------
seir <- pomp(
          data = data.frame(time = seq(1, 729, by = 7), cases = data),  #"cases" is the dataset, "time" is the observation time
          times = "time",
          t0 = 1,                      # initial time point
          dmeasure = Csnippet(dmeas),  # evaluates the density of the measurement process
          rmeasure = Csnippet(rmeas),  # simulates from the measurement process
          rprocess = SEIR_sim,          # simulates from the latent process
          statenames = c("S", "E", "I", "R"),  #state space variable name
          paramnames = c("beta", "gamma", "mu", "rho", "theta1", "theta2", "theta3"), #parameters name
          initializer = function(params, t0, ...) {
                    #Initial proportion of S, E, I, R
                    ps <- exp(params["theta1"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
                    pe <- exp(params["theta2"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
                    pi <- 1 / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
                    pr <- exp(params["theta3"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))  
                    return(setNames(as.numeric(rmultinom(1, 500, prob = c(ps, pe, pi, pr))), c("S", "E", "I", "R")))
          },
          params = c(
                    beta = log(0.00072),
                    gamma = log(0.5),
                    mu = log(0.2),
                    rho = log(1 / 4),
                    theta1 = log(90),
                    theta2 = log(2),
                    theta3 = log(7)
          )
)

## ----echo=T,eval=T-------------------------------------------------------
trans<-function(a,b,c){
          a_t<-exp(a)
          b_t<-exp(b)
          c_t<-exp(c)
          return(a_t*b_t*c_t/((1+a_t+b_t+c_t)^4))
} #Jacobian for (theta1, theta2, theta3)

seir.dprior <- function(params, ..., log) {
          f <- (dgamma(exp(params[1]), 1, 10000, log = TRUE) + params[1] +  #log prior for beta
                          dgamma(exp(params[2]), 1, 11, log = TRUE) + params[2] +  #log prior for gamma
                          dgamma(exp(params[3]), 3.2, 100, log = TRUE) + params[3] +   #log prior for mu
                          dbeta(exp(params[4]) / (1 + exp(params[4])), 3.5, 6.5, log = TRUE) +
                          params[4] - log((1 + exp(params[4])) ^ 2) +      #log prior for rho
                          log(ddirichlet(c(exp(params[5]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                                           exp(params[6]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                                           1 / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                                           exp(params[7]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7]))),
                                         c(100, 0.01, 0.4, 0.01)) * trans(params[5], params[6], params[7]))  #log prior for initial SEIR value
          )
          if (log) {
                    f
          } else {
                    exp(f)
          }
}

## ----echo=T,eval=T-------------------------------------------------------
param.initial<- c(
          beta = log(abs(rnorm(1, 0.00005, 1e-6))),
          gamma = log(abs(rnorm(1, 0.15, 0.01))),
          mu = log(abs(rnorm(1, 0.04, 0.001))),
          rho = -1,
          theta1 = 5,
          theta2 = 1,
          theta3 = 1
)

## ----echo=T,eval=T,warning=F---------------------------------------------
pmcmc1 <- pmcmc(
          pomp(seir, dprior = seir.dprior),
          #given the prior function
          start = param.initial,
          #given the initial value of the parameters
          Nmcmc = 10,
          #number of mcmc steps
          Np = 200,
          max.fail = Inf,
          proposal = mvn.rw.adaptive(0.5 * c(
                    beta = 0.1,
                    #sampling variance for beta
                    gamma = 0.1,
                    #sampling variance for gamma
                    mu = 0.1,
                    #sampling variance for mu
                    rho = 0.1,
                    #sampling variance for rho
                    theta1 = 0.1,
                    #sampling variance for theta1
                    theta2 = 0.1,
                    #sampling variance for theta2
                    theta3 = 0.1 #sampling variance for theta3
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

