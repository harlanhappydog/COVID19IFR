############################################################
#### Bayesian adjustment for preferential testing ####
#### in estimating the COVID-19 infection fatality rate ####
#### code written by Harlan Campbell ####
#### contact harlan.campbell@stat.ubc.ca ####

#### Code to reproduce illustrative example:

# Determine missing packages and load them:
required_packages <- c("quadprog" , "MASS", "xtable")
not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]    
if(length(not_installed)) install.packages(not_installed)                                           
suppressWarnings(lapply(required_packages, require, character.only = TRUE))

# functions:
invlog <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
cloglog <- function(x){log(-log(1-x))}
icloglog <- function(x){1 - exp(-exp(x))}

############
### simulate_data: function for simulating some artificial data: ###
############

simulate_data <- function(K ,
                          kprime, 
                          theta , 
                          tausquared , 
                          beta0 , 
                          sigmasquared , 
                          pop_mu , 
                          testingrate_lower , 
                          testingrate_upper , 
                          gamma_vec ){
  
  ## parameters: ##
  cloglogIFR  <- rnorm(K, mean = theta, sd=sqrt(tausquared))
  IFR         <- icloglog(cloglogIFR)
  cloglogIR   <- rnorm(K, beta0, sqrt(sigmasquared)) 
  IR          <- icloglog(cloglogIR)
  testingrate <- runif(K, testingrate_lower, testingrate_upper)
  
  ## data: ##
  kvec        <- c(1:K); population  <- rep(0,K)
  if(kprime > 0){ population[1:kprime] <- 10+rnegbin(kprime, 0.1*pop_mu-10, 1); 
  if(kprime < K){
    population[(kprime+1):K] <- 10+rnegbin(length((kprime+1):K), pop_mu-10, 1) }}
  if(kprime == 0){ population[1:K] <-  10+rnegbin(K, pop_mu-10, 1) }  
  
  tests       <- round(testingrate*population)
  cases       <- round(IR*population)
  deaths      <- rbinom(K, cases, IFR)
  
  if(kprime > 0) {
    CC_fixed <-  apply(
      cbind(1:kprime), 1, function(k){ rWNCHypergeo(1, m1=cases[k], m2=population[k]-cases[k], n=tests[k], odds = 1)})}
  
  phi_list <- confirmed_cases <- list()
  for(jj in 1:length(gamma_vec)){
    
    phi <- rep(0, K)  
    if( kprime > 0 ){ phi[1:kprime] <- 1; 
    if(kprime < K ){
      phi[(kprime+1):K] <- (seq(1, 1+gamma_vec[jj], length.out=length((kprime+1):K)))
    }}
    if( kprime == 0 ){ phi[1:K] <- 1 + seq(0, gamma_vec[jj],length.out=K) }   
    
    confirmed_cases[[jj]] <- rep(0,K);
    if( kprime > 0 ){ confirmed_cases[[jj]][1:kprime] <- CC_fixed}
    
    if(kprime < K ){
      confirmed_cases[[jj]][(kprime+1):K] <- apply(
        cbind((kprime+1):K), 1, function(k){ rWNCHypergeo(1, m1=cases[k], m2=population[k]-cases[k], n=tests[k], odds = phi[k])})
    }
    phi_list[[jj]] <- phi
  }
  
  full_data <- list(df=data.frame(kvec, population, tests, deaths, cases, IR, IFR), CC=confirmed_cases, phi= phi_list)
  
  return(full_data)}

############
### end of simulate_data function
############


#########
set.seed(1234)
OD <- simulate_data(K = 12, kprime=0, theta = cloglog(0.02),  tausquared = 0.005, beta0 = cloglog(0.20), sigmasquared = 0.25, pop_mu = 2000, testingrate_lower = 0.05, testingrate_upper = 0.10, gamma_vec=c(0,4,11,22))

print(cbind(OD$df[,c(2:4)], CC=matrix(unlist(OD$CC),12,), phi=matrix(unlist(OD$phi),12,), OD$df[,c(5,6,7)]), digits=3)
xtable(cbind(OD$df[,c(2:4)], CC=matrix(unlist(OD$CC),12,), phi=matrix(unlist(OD$phi),12,), OD$df[,c(5,6,7)]), digits=3)
