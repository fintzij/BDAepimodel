## ----echo=T,eval=T-------------------------------------------------------
library(pomp)
library(ggplot2)
library(gtools)
set.seed(52787)
st.date <- as.Date("1978/1/22")
en.date <- as.Date("1978/2/4")
date.seq <- seq(st.date,en.date, by = "day")

bbs <- data.frame(time = rep(date.seq,1),
                  count = c(rep("Observed count", 14)),
                  I = c(3,8,28,76,222,293,257,237,192,126,70,28,12,5))

## ----echo=T,eval=T-------------------------------------------------------
rmeas <- "
cases=rnbinom_mu(exp(phi),exp(rho)*I/(1+exp(rho)));    //represent the data
"
dmeas<-"
lik=dnbinom_mu(cases,exp(phi),exp(rho)*I/(1+exp(rho)),give_log); //return the loglikelihood
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
" 

## ----echo=T,eval=T-------------------------------------------------------
sir <- pomp(
  data = data.frame(cases = bbs[,3], time = seq(0, 13, by = 1)),  #"cases" is the dataset, "time" is the observation time
  times = "time",
  t0 = 0,  #initial time point
  dmeasure = Csnippet(dmeas),  
  rmeasure = Csnippet(rmeas),  #return the likelihood
  rprocess = euler.sim(step.fun = Csnippet(sir.step), delta.t = 1/12),
  statenames = c("S", "I", "R"),  #state space variable name
  paramnames = c("bet", "mu", "rho", "theta1", "theta2", "phi"), #parameters name
  initializer = function(params, t0, ...) {
    ps <-exp(params["theta1"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  #given the initial proportion of susceptible
    pi <- 1 / (1 + exp(params["theta1"]) + exp(params["theta2"]))                     #given the initial proportion of infection
    pr <-exp(params["theta2"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]))  #given the initial proportion of recovery
    return(setNames(as.numeric(rmultinom(
      1, 763, prob = c(ps, pi, pr)
    )), c("S", "I", "R")))
  },
  params = c(    #for data generation(We didn't generate data here, but pomp need this part.)
    bet = -5,  
    mu = -1,
    rho = 3,
    theta1 = 5,
    theta2 = 1,
    phi = 2
  )
)

## ----echo=T,eval=T-------------------------------------------------------
sir.dprior <- function(params, ..., log) {
  f <-
    (
      dgamma(
        exp(params[1]),
        shape = 0.001,
        rate = 1,
        log = TRUE
      ) + params[1] +          #log prior for log(infection rate) "beta"
        dgamma(
          exp(params[2]),
          shape = 1,
          rate = 2,
          log = TRUE
        ) + params[2] +        #log prior for log(recovery rate) "mu"
        dbeta(exp(params[3]) / (1 + exp(params[3])), 2, 1, log = TRUE) +
        params[3] - log((1 + exp(params[3])) ^ 2) + #log prior for logit(sampling probablity) "rho"
        log(ddirichlet(c(
          exp(params[4]) / (1 + exp(params[4]) + exp(params[5])),
          1 / (1 + exp(params[4]) + exp(params[5])),
          exp(params[5]) / (1 + exp(params[4]) + exp(params[5]))
        ), c(900, 3, 9))) + params[4] + params[5] - 3 * log(1 + exp(params[4]) +
                                                              exp(params[5])) + #log prior for logit(initial value)
        dgamma(exp(params[6]), 1, 0.1, log = TRUE) + params[6]  #log prior for log(negative size parameter) "phi"
    )
  if (log)
    f
  else
    exp(f)
} 

## ----echo=T,eval=T-------------------------------------------------------
param.initial <- c(  #initial value of the parameters of the PMCMC
  bet = -6,
  mu = -1,
  rho = 2.5,
  theta1 = 6,
  theta2 = 1,
  phi = 2
)

## ----echo=T,eval=T,warning=F---------------------------------------------
pmcmc1 <- pmcmc(
          pomp(sir, dprior = sir.dprior),
          #given the prior function
          start = param.initial,
          #given the initial value of the parameters
          Nmcmc = 10,
          #number of mcmc steps
          Np = 500,
          max.fail = Inf,
          proposal = mvn.rw.adaptive(0.1 * c(
                    bet = 0.1,
                    #sampling variance for beta
                    mu = 0.1,
                    #sampling variance for mu
                    rho = 0.3,
                    #sampling variance for rho
                    theta1 = 0.2,
                    #sampling variance for theta1
                    theta2 = 0.2,
                    #sampling variance for theta2
                    phi = 0.5   #sampling variance for phi
          ),
          shape.start = 100, 
          scale.start = 100,
          max.scaling = 2)
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

end_time <- Sys.time();
run_time <- difftime(end_time, start_time, units = "hours")

pomp_results <- list(time = run_time, results = pmcmc1) 

## ----echo=T,eval=T-------------------------------------------------------
rmeas <- "
cases=rnbinom_mu(exp(phi),exp(rho)*I/(1+exp(rho)));    //represent the data
"
dmeas<-"
lik=dnbinom_mu(cases,exp(phi),exp(rho)*I/(1+exp(rho)),give_log); //return the loglikelihood
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

## ----echo=T,eval=T-------------------------------------------------------
seir <- pomp(
          data = data.frame(cases = bbs[,3], time = seq(0, 13, by = 1)),  #"cases" is the dataset, "time" is the observation time
          times = "time",
          t0 = 0,  #initial time point
          dmeasure = Csnippet(dmeas),  
          rmeasure = Csnippet(rmeas),  #return the likelihood
          rprocess = euler.sim(step.fun = Csnippet(seir.step), delta.t = 1/12), #delta.t is here
          statenames = c("S", "E", "I", "R"), #state space variable name
          paramnames = c("beta", "gamma", "mu", "rho", "theta1", "theta2", "theta3", "phi"), #parameters name
          initializer = function(params, t0, ...) {
                    #Initial proportion of S, E, I, R
                    ps <- exp(params["theta1"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
                    pe <- exp(params["theta2"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
                    pi <- 1 / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))
                    pr <- exp(params["theta3"]) / (1 + exp(params["theta1"]) + exp(params["theta2"]) + exp(params["theta3"]))  
                    return(setNames(as.numeric(rmultinom(1, 763, prob = c(ps, pe, pi, pr))), c("S", "E", "I", "R")))
          },
          params = c(    #for data generation(We didn't generate data here, but pomp need this part.)
                    beta = -6,
                    gamma = 1.5,
                    mu = -1,
                    rho = 2.5,
                    theta1 = 6,
                    theta2 = 1,
                    theta3 = 1,
                    phi = 2
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
          f <- (dgamma(exp(params[1]), 0.001, 1, log = TRUE) + params[1] +  #log prior for beta
                          dgamma(exp(params[2]), 0.001, 1, log = TRUE) + params[2] +  #log prior for gamma
                          dgamma(exp(params[3]), 1, 2, log = TRUE) + params[3] +   #log prior for mu
                          dbeta(exp(params[4]) / (1 + exp(params[4])), 2, 1, log = TRUE) +
                          params[4] - log((1 + exp(params[4])) ^ 2) +      #log prior for rho
                          log(ddirichlet(c(exp(params[5]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                                           exp(params[6]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                                           1 / (1 + exp(params[5]) + exp(params[6]) + exp(params[7])),
                                           exp(params[7]) / (1 + exp(params[5]) + exp(params[6]) + exp(params[7]))),
                                         c(900, 6, 3, 9)) * trans(params[5], params[6], params[7])) + # log prior for p_t1
                          dgamma(exp(params[6]), 1, 0.1, log = TRUE) + params[6]  #log prior for log(negative size parameter) "phi"
          )
          if (log) {
                    f
          } else {
                    exp(f)
          }
}

## ----echo=T,eval=T-------------------------------------------------------
param.initial <- c(  #initial value of the parameters of the PMCMC
  beta = -6,
  gamma = 2,
  mu = -1,
  rho = 2.5,
  theta1 = 6,
  theta2 = 1,
  theta3 = 1,
  phi = 2
)

## ----echo=T,eval=T,warning=F---------------------------------------------
pmcmc1 <- pmcmc(
          pomp(seir, dprior = seir.dprior),
          #given the prior function
          start = param.initial,
          #given the initial value of the parameters
          Nmcmc = 2000,
          #number of mcmc steps
          Np = 10,
          max.fail = Inf,
          proposal = mvn.rw.adaptive(0.1 * c(
                    beta = 0.1, #sampling variance for beta
                    gamma = 0.1, # proposal variance for gamma
                    mu = 0.1, #sampling variance for mu
                    rho = 0.3, #sampling variance for rho
                    theta1 = 0.2, #sampling variance for theta1
                    theta2 = 0.2, #sampling variance for theta2
                    theta3 = 0.2,
                    phi = 0.5   #sampling variance for phi
          ),
          shape.start = 100, 
          scale.start = 100,
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

end_time <- Sys.time();
run_time <- difftime(end_time, start_time, units = "hours")

pomp_results <- list(time = run_time, results = pmcmc1) 

