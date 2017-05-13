## ----echo=T,eval=T-------------------------------------------------------
library(BDAepimodel)
library(coda)
library(gtools)
library(Rcpp)
library(pomp)
library(ggplot2)

## ----echo=T,eval=T-------------------------------------------------------
data<-c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,2,0,0,0,1,0,1,0,0,0,0,0,2,2,0,3,0,1,0,0,0,0,0,2,1,0,0,0,0,0,0,0,1,0,0,0,3,1,0,0,1,2,2,0,4,0,3,2,0,2,5,2,1,2,2,2,1,1,3,1,2,5,1,2,1,0,1,1,4,1,0,1,0,1,2,0,2,0,1,0,1,0,3,1,1,2,1,1,1,2,3,0,1,2,1,3,3,2,2,0,0,3,2,3,3,0,2,0,2,1,2,3,2,1,0,1,1,0,3,1,2,2,1,1,1,1,0,1,0,0,1,0,1,1,1,1,0,0,0,1,0,0,1,1,1,0,0,0,0,1,0,0,0,2,0,0,0,1,0,0,0,0,0,1,0,3,0,1,1,0,0,0,1,1,1,1,0,0,2,0,1,1,1,0,1,0,0,1,1,0,1,0,1,1)

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
rate[1]=exp(gamma);     //Rate of S to E
rate[2]=exp(mu);      //Rate of E to I
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
  data = data.frame(cases = data, time = seq(1, 1457, by = 7)), #"cases" is the dataset, "time" is the observation time
  times = "time",
  t0 = 1,  #initial time point
  dmeasure = Csnippet(dmeas),
  rmeasure = Csnippet(rmeas),
  rprocess = euler.sim(step.fun = Csnippet(seir.step), delta.t = 1/4), #delta.t is here
  statenames = c("S", "E", "I", "R"), #state space variable name
  paramnames = c("beta", "gamma", "mu", "rho", "theta1", "theta2", "theta3"), #parameters name
  initializer = function(params, t0, ...) {
            ps <-exp(params["theta1"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
            pe <- exp(params["theta2"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
            pi <- 1 / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
            pr <- exp(params["theta3"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))  #Initial proportion of S, E, I, R
    return(setNames(as.numeric(rmultinom(1, 400, prob = c(ps, pe, pi, pr))), c("S", "E", "I", "R")))
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

## ----echo=T,eval=T,warning=F---------------------------------------------
trans<-function(a,b,c){
  a_t<-exp(a)
  b_t<-exp(b)
  c_t<-exp(c)
  return(a_t*b_t*c_t/((1+a_t+b_t+c_t)^4))
} #Jacobian for (theta1, theta2, theta3)

seir.dprior <- function(params, ..., log) {
  f <- (
    dgamma(exp(params[1]), 0.6, 10000, log = TRUE) + params[1] +  #log prior for beta
      dgamma(exp(params[2]), 1.5, 25, log = TRUE) + params[2] +  #log prior for gamma
      dgamma(exp(params[3]), 1.25, 70, log = TRUE) + params[3] +   #log prior for mu
      dbeta(exp(params[4]) / (1 + exp(params[4])), 5, 50, log = TRUE) +
      params[4] - log((1 + exp(params[4])) ^ 2) +      #log prior for rho
      log(ddirichlet(c(exp(params[5]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                       exp(params[6]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                       1 / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                       exp(params[7]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7]))),
                     c(90, 0.5, 0.5, 0.01)) * trans(params[5], params[6], params[7]))  #log prior for initial SEIR value
  )
  if (log) {
          f  
  } else {
          exp(f)         
  }
}

## ----echo=T,eval=T,warning=F---------------------------------------------
param.initial<- c(
  beta = log(abs(rnorm(1, 0.000055, 1e-6))),
  gamma = log(abs(rnorm(1, 0.03, 1e-3))),
  mu = log(abs(rnorm(1, 0.01, 0.01))),
  rho = -2,
  theta1 = 6,
  theta2 = -1,
  theta3 = -2
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
                    0.1 * c(
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
                    scale.start = 1000,
                    shape.start = 1000,
                    max.scaling = 1.2
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
# measurement process
rmeas<-"
cases=rbinom(I,exp(rho)/(1+exp(rho)));//represent the data
"
dmeas<-"
if(cases > 0 && I == 0) {
          if(give_log) {
                    lik = R_NegInf;
          } else {
                    lik = 0;
          }
} else {
          lik=dbinom(cases,I,exp(rho)/(1+exp(rho)),give_log); //return the loglikelihood
}
"

## ----echo=T,eval=T-------------------------------------------------------
sir.step<-"
double rate[2];
double dN[2];
rate[0]=exp(bet)*I;    //Infection rate
rate[1]=exp(mu);       //recovery rate
reulermultinom(1,S,&rate[0],dt,&dN[0]);   //generate the number of newly infected people
reulermultinom(1,I,&rate[1],dt,&dN[1]);   //generate the number of newly recovered people
if(!R_FINITE(S)) Rprintf(\"%lg %lg %lg %lg %lg %lg %lg %lg %lg\\n\",dN[0],rate[0],dN[1],rate[1],bet,mu,S,I,R);
S+=-dN[0];          //update the number of Susceptible
I+=dN[0]-dN[1];     //update the number of Infection
R+=dN[1];           //update the number of Recovery
if(I<0) Rprintf(\"%lg %lg %lg %lg %lg %lg %lg %lg %lg\\n\",dN[0],rate[0],dN[1],rate[1],bet,mu,S,I,R);
" 

## ----echo=T,eval=T-------------------------------------------------------
sir <- pomp(
          data = data.frame(cases = data, time = seq(1, 1457, by = 7)),  #"cases" is the dataset, "time" is the observation time
          times = "time",
          t0 = 1,  #initial time point
          dmeasure = Csnippet(dmeas),  
          rmeasure = Csnippet(rmeas),  #return the likelihood
          rprocess = euler.sim(step.fun = Csnippet(sir.step), delta.t = 1/3), # delta.t is here
          statenames = c("S", "I", "R"),  #state space variable name
          paramnames = c("bet", "mu", "rho", "theta1", "theta2"), #parameters name
          initializer = function(params, t0, ...) {
                    ps <-exp(params["theta1"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  # initial prob of susceptible
                    pi <- 1 / (1 + exp(params["theta1"]) + exp(params["theta2"]))                     # initial prob of infection
                    pr <-exp(params["theta2"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  # initial prob of recovery
                    return(setNames(as.numeric(rmultinom(
                              1, 400, prob = c(ps, pi, pr)
                    )), c("S", "I", "R")))
          },
          params = c(    #for data generation(We didn't generate data here, but pomp need this part.)
                    bet = -5,  
                    mu = -2,
                    rho = 0,
                    theta1 = 6,
                    theta2 = 1
          )
)

## ----echo=T,eval=T-------------------------------------------------------
sir.dprior <- function(params, ..., log) {
          f <- (dgamma(exp(params[1]), shape = 0.6, rate = 10000, log = TRUE) + params[1] +#log prior for log(infection rate) "beta"
                          dgamma(exp(params[2]), shape = 1.25, rate = 70, log = TRUE) + 
                          params[2] + #log prior for log(recovery rate) "mu"
                          
                          dbeta(exp(params[3]) / (1 + exp(params[3])), 5, 50, log = TRUE) +
                          params[3] - log((1 + exp(params[3])) ^ 2) + #log prior for logit(sampling probablity) "rho"
                          
                          log(ddirichlet(c(exp(params[4]) / (1 + exp(params[4]) + exp(params[5])),
                                           1 / (1 + exp(params[4]) + exp(params[5])),
                                           exp(params[5]) / (1 + exp(params[4]) + exp(params[5]))
                          ), c(90, 0.5, 0.01))) +
                          params[4] + params[5] - 3 * log(1 + exp(params[4]) + exp(params[5])) #log prior of logit(initial value)
          )
          if (log) {
                    f
          } else {
                    exp(f)
          }
} 



## ----echo=T,eval=T-------------------------------------------------------
param.initial <- c(
          bet = log(abs(rnorm(1, 0.00007, 1e-5))),
          mu = log(abs(rnorm(1, 0.005, 0.0001))),
          rho = -2,
          theta1 = 6,
          theta2 = -1
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
                    bet = 0.1,
                    #sampling variance for beta
                    mu = 0.1,
                    #sampling variance for mu
                    rho = 0.1,
                    #sampling variance for rho
                    theta1 = 0.1,
                    #sampling variance for theta1
                    theta2 = 0.1 #sampling variance for theta2
          ), 
          scale.start = 100, 
          shape.start = 100,
          max.scaling = 1.2)
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

