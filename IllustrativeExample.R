############################################################
#### Bayesian adjustment for preferential testing ####
#### in estimating the COVID-19 infection fatality rate ####
#### code written by Harlan Campbell ####
#### contact harlan.campbell@stat.ubc.ca ####

#### Code to reproduce illustrative example analysis:

# Determine missing packages and load them:
required_packages <- c("xtable", "MCMCvis", "rjags", "coda", "MCMCpack", "BiasedUrn", "rvest", "runjags", "MASS", "HDInterval", "gridGraphics", "gridExtra", "forestplot", "latex2exp")
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


############
### JAGS Model for "M_{gamma=0}" ###
############
model<- '
model{	
#likelihood
for (k in 1:K){
cloglog_IFR[k] ~ dnorm(theta, inv.var_tau);
cloglog_infectionrate[k] ~ dnorm(beta0, inv.var_sig) ;
cloglog(IFR[k]) <- cloglog_IFR[k];
cloglog(infectionrate[k]) <- cloglog_infectionrate[k];

confirmed_cases[k] ~ dhyper(cases[k], population[k]-cases[k], tests[k], phi[k]);
cases[k] ~ dbin(infectionrate[k], population[k]);
deaths[k] ~ dbin(IFR[k], cases[k]);
}

#priors
for (k in 1:K){ phi[k]<-1;}

icloglog_theta ~ dunif(0, 1); 
icloglog_beta0 ~ dunif(0, 1);
theta <- log(-log(1-icloglog_theta));
beta0 <- log(-log(1-icloglog_beta0));

inv.var_sig   <- (1/sd_sig)^2 ;
sd_sig     ~ dnorm(0, 1/1) T(0,);
inv.var_tau   <- (1/sd_tau)^2 ;
sd_tau     ~ dnorm(0, 1/0.01) T(0,);
}'
#    
cat(model, file="JAGS_ignore.txt")
############
### end of JAGS Model for "M_{gamma=0}" ###
############

############
### JAGS model for "M_{1}" ###
############
model<- '
model{	
#likelihood
for (k in 1:K){
cloglog_IFR[k] ~ dnorm(theta, inv.var_tau);
cloglog_infectionrate[k] ~ dnorm(beta0, inv.var_sig) ;
cloglog(IFR[k]) <- cloglog_IFR[k];
cloglog(infectionrate[k]) <- cloglog_infectionrate[k];

confirmed_cases[k] ~ dhyper(cases[k], population[k]-cases[k], tests[k], phi[k]);
cases[k] ~ dbin(infectionrate[k], population[k]);
deaths[k] ~ dbin(IFR[k], cases[k]);
}

#priors
for (k in 1:K){ phi[k] ~ dunif(1, 1+gamma);}
gamma  ~ dexp(0.5);

icloglog_theta ~ dunif(0, 1); 
icloglog_beta0 ~ dunif(0, 1);
theta <- log(-log(1-icloglog_theta));
beta0 <- log(-log(1-icloglog_beta0));

inv.var_sig   <- (1/sd_sig)^2 ;
sd_sig     ~ dnorm(0, 1/1) T(0,);
inv.var_tau   <- (1/sd_tau)^2 ;
sd_tau     ~ dnorm(0, 1/0.01) T(0,);
}'
#    
cat(model, file="JAGS_unif.txt")
############
### END of JAGS model for "M_{1}" ###
############


############
### illustrative: function for running M_{1} and M_{gamma=0} models for all values of gamma ###
############

illustrative <- function(observed_data, MCMCiter=50000){
  
  ###
  datalist <- list(K=length(observed_data$kvec), confirmed_cases=observed_data$confirmed_cases, deaths=observed_data$deaths, population=observed_data$population, tests=observed_data$tests)
  
  jags.m_unif <- jags.model(file = "JAGS_unif.txt", data = datalist, n.chains = 5, n.adapt = 50000, inits= NULL)
  
  ###
  datalist <- list(K=length(observed_data$kvec), confirmed_cases=observed_data$confirmed_cases, deaths=observed_data$deaths, population=observed_data$population, tests=observed_data$tests)
  
  jags.m_ignore <- jags.model(file = "JAGS_ignore.txt", data = datalist, n.chains = 5, n.adapt = 50000, inits= NULL)
  
  ######
  params <- c("phi", "theta", "beta0",  "IFR", "infectionrate", "sd_sig",  "gamma", "sd_tau")		
  
  samps_unif <- coda.samples(jags.m_unif, params, n.iter=MCMCiter,  n.burnin = MCMCiter*0.2, thin=50)
  
  samps_ignore <- coda.samples(jags.m_ignore, params[-which.max(params=="gamma")], n.iter=MCMCiter,  n.burnin = MCMCiter*0.2, thin=50)
  
  median_unif <- summary(samps_unif)$quantiles[,c(3)]
  median_ignore <- summary(samps_ignore)$quantiles[,c(3)]
  
  HDI<-hdi(samps_unif)
  QQ_unif<-t(rbind(HDI[1,],median_unif,HDI[2,]))
  HDI<-hdi(samps_ignore)
  QQ_ignore<-t(rbind(HDI[1,],median_ignore,HDI[2,]))
  
  QQA<-cbind(QQ_unif, QQ_ignore[match(rownames(QQ_unif), rownames(QQ_ignore)),])
  goodQ<-QQA[c("theta", "beta0",   "sd_tau", "sd_sig" , "gamma"),]
  
  return(list(goodQ = goodQ, QQA=QQA, samps_unif = samps_unif, samps_ignore = samps_ignore))
}
############
### end of illustrative function ###
############



### START OF CREATING DATA AND RESULTS ###

#########
set.seed(1234)
OD <- simulate_data(K = 12, kprime=0, theta = cloglog(0.02),  tausquared = 0.005, beta0 = cloglog(0.20), sigmasquared = 0.25, pop_mu = 2000, testingrate_lower = 0.05, testingrate_upper = 0.10, gamma_vec=c(0,4,11,22))

print(cbind(OD$df[,c(2:4)], CC=matrix(unlist(OD$CC),12,), phi=matrix(unlist(OD$phi),12,), OD$df[,c(5,6,7)]), digits=3)

xtable(cbind(OD$df[,c(2:4)], CC=matrix(unlist(OD$CC),12,),  OD$df[,c(5,6,7)], phi=round(matrix(unlist(OD$phi),12,)[,-1],2)), digits=3)

### "When $\gamma=0$, the number of true cases (i.e. actual infections) is approximately 14 times higher than the number of confirmed cases.  In contrast, when $\gamma=22$, the number of true cases is only about 5 times higher than the number of confirmed cases." ###

mean(OD$df[,"cases"]/OD$CC[[1]])
mean(OD$df[,"cases"]/OD$CC[[4]])

##############################

### "Each model is fit using JAGS (just another Gibbs' sampler) \citep{kruschke2014doing}, with 5 independent chains, each with 500,000 draws (20\% burnin, thinning of 50)." ###  

MCMCiter_sim <- 500000


####### results for gamma = 0
results0  <- illustrative(cbind(OD$df, confirmed_cases =unlist(OD$CC[[1]])),  MCMCiter = MCMCiter_sim)
print(results0$goodQ,digits=2)

MCMCtrace(results0$samps_unif, params= c("theta","beta0", "gamma"), priors=cbind(cloglog(runif(MCMCiter_sim)),cloglog(runif(MCMCiter_sim)),rexp(MCMCiter_sim,0.5)) , main_den = c(			   
  TeX("Density $\\theta$"),
  TeX("Density $\\beta$"),
  TeX("Density $\\gamma$")),
  main_tr = c(					   					                         
    TeX("Trace $\\theta$"),
    TeX("Trace $\\beta$"),
    TeX("Trace $\\gamma$")),
  filename= "MCMC_theta_gamma0.pdf")


MCMCtrace(results0$samps_unif, ISB = FALSE, params= c("phi\\[1\\]","phi\\[6\\]","phi\\[12\\]"), priors=apply(cbind(rexp(MCMCiter_sim,0.5)),1,function(x) runif(1, 1, 1+x) ) , main_den = c(			   
  TeX("Density $\\phi_{1}$"),
  TeX("Density $\\phi_{6}$"),
  TeX("Density $\\phi_{12}$")
), 
main_tr = c(					   					                         
  TeX("Trace $\\phi_{1}$"),
  TeX("Trace $\\phi_{6}$"),
  TeX("Trace $\\phi_{12}$")), 
filename= "MCMC_phi_gamma0.pdf")


####### results for gamma = 4

results4  <- illustrative(cbind(OD$df, confirmed_cases =unlist(OD$CC[[2]])),  MCMCiter = MCMCiter_sim)
print(results4$goodQ,digits=2)

MCMCtrace(results4$samps_unif, params= c("theta","beta0", "gamma"), priors=cbind(cloglog(runif(MCMCiter_sim)),cloglog(runif(MCMCiter_sim)),rexp(MCMCiter_sim,0.5)) , main_den = c(		   
  TeX("Density $\\theta$"),
  TeX("Density $\\beta$"),
  TeX("Density $\\gamma$")),
  main_tr = c(					   					                         
    TeX("Trace $\\theta$"),
    TeX("Trace $\\beta$"),
    TeX("Trace $\\gamma$")),
  filename= "MCMC_theta_gamma4.pdf")

MCMCtrace(results4$samps_unif, ISB = FALSE, params= c("phi\\[1\\]","phi\\[6\\]","phi\\[12\\]"), priors=apply(cbind(rexp(MCMCiter_sim,0.5)),1,function(x) runif(1, 1, 1+x) ) , main_den = c(			   
  TeX("Density $\\phi_{1}$"),
  TeX("Density $\\phi_{6}$"),
  TeX("Density $\\phi_{12}$")), 
  main_tr = c(					   					                         
    TeX("Trace $\\phi_{1}$"),
    TeX("Trace $\\phi_{6}$"),
    TeX("Trace $\\phi_{12}$")), 
  filename= "MCMC_phi_gamma4.pdf")

####### results for gamma = 11

results11  <- illustrative(cbind(OD$df, confirmed_cases =unlist(OD$CC[[3]])),  MCMCiter = MCMCiter_sim)
print(results11$goodQ,digits=2)

MCMCtrace(results11$samps_unif, params= c("theta","beta0", "gamma"), priors=cbind(cloglog(runif(MCMCiter_sim)),cloglog(runif(MCMCiter_sim)),rexp(MCMCiter_sim,0.5)) , main_den = c(		   
  TeX("Density $\\theta$"),
  TeX("Density $\\beta$"),
  TeX("Density $\\gamma$")),
  main_tr = c(					   					                         
    TeX("Trace $\\theta$"),
    TeX("Trace $\\beta$"),
    TeX("Trace $\\gamma$")),
  filename= "MCMC_theta_gamma11.pdf")

MCMCtrace(results11$samps_unif, ISB = FALSE, params= c("phi\\[1\\]","phi\\[6\\]","phi\\[12\\]"), priors=apply(cbind(rexp(MCMCiter_sim,0.5)),1,function(x) runif(1, 1, 1+x) ) , main_den = c(			   
  TeX("Density $\\phi_{1}$"),
  TeX("Density $\\phi_{6}$"),
  TeX("Density $\\phi_{12}$")), 
  main_tr = c(					   					                         
    TeX("Trace $\\phi_{1}$"),
    TeX("Trace $\\phi_{6}$"),
    TeX("Trace $\\phi_{12}$")), 
  filename= "MCMC_phi_gamma11.pdf")

####### results for gamma = 22

results22  <- illustrative(cbind(OD$df, confirmed_cases =unlist(OD$CC[[4]])),  MCMCiter = MCMCiter_sim)
print(results22$goodQ,digits=2)

MCMCtrace(results22$samps_unif, params= c("theta","beta0", "gamma"), priors=cbind(cloglog(runif(MCMCiter_sim)),cloglog(runif(MCMCiter_sim)),rexp(MCMCiter_sim,0.5)) , main_den = c(		
  TeX("Density $\\theta$"),
  TeX("Density $\\beta$"),
  TeX("Density $\\gamma$")),
  main_tr = c(					   					                         
    TeX("Trace $\\theta$"),
    TeX("Trace $\\beta$"),
    TeX("Trace $\\gamma$")),
  filename= "MCMC_theta_gamma22.pdf")

MCMCtrace(results22$samps_unif, ISB = FALSE, params= c("phi\\[1\\]","phi\\[6\\]","phi\\[12\\]"), priors=apply(cbind(rexp(MCMCiter_sim,0.5)),1,function(x) runif(1, 1, 1+x) ) , main_den = c(			   
  TeX("Density $\\phi_{1}$"),
  TeX("Density $\\phi_{6}$"),
  TeX("Density $\\phi_{12}$")
), 
main_tr = c(					   					                         
  TeX("Trace $\\phi_{1}$"),
  TeX("Trace $\\phi_{6}$"),
  TeX("Trace $\\phi_{12}$")), 
filename= "MCMC_phi_gamma22.pdf")


#### Saving results

results_list<-list(OD, results0, results4, results11, results22)
saveRDS(results_list, "illus_example.rds")

results_list<-readRDS("illus_example.rds")

results0<-results_list[[2]]
results4<-results_list[[3]]
results11<-results_list[[4]]
results22<-results_list[[5]]
##############################
summ <- rbind(
  results0$QQA[c("theta", "beta0", "sd_tau", "sd_sig", "gamma"),],
  results4$QQA[c("theta", "beta0", "sd_tau", "sd_sig", "gamma"),],
  results11$QQA[c("theta", "beta0", "sd_tau", "sd_sig", "gamma"),],
  results22$QQA[c("theta", "beta0", "sd_tau", "sd_sig", "gamma"),])

summ0 <- summ

summ0[rownames(summ0)%in%"theta",] <- icloglog(summ[rownames(summ)%in%"theta",])
summ0[rownames(summ0)%in%"beta0",] <- icloglog(summ[rownames(summ)%in%"beta0",])

print(round(summ0,2))
xtable(summ0[], digits=3)

