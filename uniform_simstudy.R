####################################################################################
# Determine missing packages and load them:
required_packages <- c("nimble", "MASS", "BiasedUrn", "xtable", "HDInterval")
not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]    
if(length(not_installed)) install.packages(not_installed)                                           
suppressWarnings(lapply(required_packages, require, character.only = TRUE))

############
invlog <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}
cloglog <- function(x){log(-log(1-x))}
icloglog <- function(x){1 - exp(-exp(x))}
lseq <- function(from, to, length.out){10^(seq(log10(from), log10(to), length.out = length.out))}


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
        phi[(kprime+1):K] <- runif(length((kprime+1):K), 1, 1+gamma_vec[jj])
#         phi[(kprime+1):K] <- seq(1, 1+gamma_vec[jj], length.out=length((kprime+1):K))
    }}
    if( kprime == 0 ){ phi[1:K] <- runif(K, 1, 1+gamma_vec[jj]) }   
#     if( kprime == 0 ){ phi[1:K] <- seq(1, 1+gamma_vec[jj], length.out=K) }    
    
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

K_sim <- 30
kprime_sim <- 5
source('initialize_MODELS.R', chdir = TRUE)
#################################################

sim_once<-function(the_K = K_sim,
                   the_kprime = kprime_sim,
                   true_pop_mu = 1000, 
                   the_gamma_vec = c(0,5,20,50),
                   the_lambda_vec = c(0.5, 1),
                   MCMCiter = 10000){
  
  OB_data <- simulate_data(K = the_K,
                           kprime = the_kprime,
                           theta = cloglog(0.02), 
                           tausquared = 0.005, 
                           beta0 = cloglog(0.2), 
                           sigmasquared = 0.25, 
                           pop_mu = true_pop_mu, 
                           testingrate_lower = 0.01, 
                           testingrate_upper = 0.10, 
                           gamma_vec = the_gamma_vec)
  HQD <- HQS <- HQ<-list()
  
  CmodelS$deaths     <- OB_data$df$deaths
  CmodelS$population <- OB_data$df$population
  CmodelS$tests      <- OB_data$df$tests
  CmodelS$confirmed_cases <- OB_data$CC[[1]]
  bb<-1
  while(length(bb)<=1) {  # in case intial values allow for good mixing
    # Reasonable initial values:
    cloglog_infectionrate_init <- runif(the_K, 0, 0.5)
    cloglog_IFR_init <- runif(the_K, 0, 0.1)
    inits <- list(cloglog_IFR = cloglog_IFR_init,
                  cloglog_infectionrate = cloglog_infectionrate_init,
                  icloglog_theta = mean(cloglog_IFR_init),
                  icloglog_beta0 = mean(cloglog_infectionrate_init),
                  sd_tau = abs(rnorm(1, sd=0.1)),
                  sd_sig = abs(rnorm(1, sd=1)))
    
    
    CmodelS$setInits(inits)
    
    mcmcsamplesS <- suppressMessages(runMCMC(CMCMC2S, thin = 50, nburnin = MCMCiter*0.2, niter = MCMCiter, progressBar = FALSE))
    bb <- unique(mcmcsamplesS[ , 'theta']);
    if(length(bb)<=1){
      print("intial values failed")
      par(mfrow=c(3,1))
      plot(mcmcsamplesS[,"sd_tau"], main= "intial values failed", type='l')
      plot(mcmcsamplesS[,"beta0"], main= "intial values failed", type='l')
      plot(mcmcsamplesS[,"theta"], main= "intial values failed", type='l')
    }
  }
  
  
  
  for(hh in 1:length(the_lambda_vec)){
    #########  M_[1:K] model with all values of gamma ###########
    
    Cmodel$deaths     <- OB_data$df$deaths
    Cmodel$population <- OB_data$df$population
    Cmodel$tests      <- OB_data$df$tests
    Cmodel$lambda     <- the_lambda_vec[hh]
    

    CmodelD$deaths     <- OB_data$df$deaths
    CmodelD$population <- OB_data$df$population
    CmodelD$tests      <- OB_data$df$tests
    CmodelD$lambda     <- the_lambda_vec[hh]    
    
     
    QD <- QS <- Q <- samplesD <- samplesS <- samples <- list();
    for(xx in 1:length(the_gamma_vec)){

      Cmodel$confirmed_cases <- OB_data$CC[[xx]]
      aa<-1
      while(length(aa)<=1) {  # in case intial values allow for good mixing
        # Reasonable initial values:
        cloglog_infectionrate_init <- runif(the_K, 0, 0.5)
        cloglog_IFR_init <- runif(the_K, 0, 0.1)
        gamma_intial <- rexp(1, the_lambda_vec[hh])
        inits <- list(cloglog_IFR = cloglog_IFR_init,
                      cloglog_infectionrate = cloglog_infectionrate_init,
                      gamma = gamma_intial,
                      phi = c(rep(1,the_kprime), runif(the_K-the_kprime, 1, 1+gamma_intial)),
                      icloglog_theta = mean(cloglog_IFR_init),
                      icloglog_beta0 = mean(cloglog_infectionrate_init),
                      sd_tau = abs(rnorm(1, sd=0.1)),
                      sd_sig = abs(rnorm(1, sd=1)))
        
        
        Cmodel$setInits(inits)
        
        mcmcsamples <- suppressMessages(runMCMC(CMCMC2, thin = 50, nburnin = MCMCiter*0.2, niter = MCMCiter, progressBar = FALSE))
        aa <- unique(mcmcsamples[ , 'gamma']);
        if(length(aa)<=1){
          print("intial values failed")
          par(mfrow=c(3,1))
          plot(mcmcsamples[,"sd_tau"], main= "intial values failed", type='l')
          plot(mcmcsamples[,"gamma"], main= "intial values failed", type='l')
          plot(mcmcsamples[,"theta"], main= "intial values failed", type='l')
        }
      }

      CmodelD$confirmed_cases <- OB_data$CC[[xx]]      
      dd<-1
      while(length(dd)<=1) {  # in case intial values allow for good mixing
        # Reasonable initial values:
        cloglog_infectionrate_init <- runif(the_K, 0, 0.5)
        cloglog_IFR_init <- runif(the_K, 0, 0.1)
        gamma_intial <- rexp(1, the_lambda_vec[hh])
        inits <- list(cloglog_IFR = cloglog_IFR_init,
                      cloglog_infectionrate = cloglog_infectionrate_init,
                      gamma = gamma_intial,
                      phi = c(rep(-99,the_kprime), runif(the_K-the_kprime, 1, 1+gamma_intial)),
                      icloglog_theta = mean(cloglog_IFR_init),
                      icloglog_beta0 = mean(cloglog_infectionrate_init),
                      sd_tau = abs(rnorm(1, sd=0.1)),
                      sd_sig = abs(rnorm(1, sd=1)))
        
        
        CmodelD$setInits(inits)
        
        mcmcsamplesD <- suppressMessages(runMCMC(CMCMC2D, thin = 50, nburnin = MCMCiter*0.2, niter = MCMCiter, progressBar = FALSE))
        dd <- unique(mcmcsamplesD[ , 'gamma']);
        if(length(aa)<=1){
          print("intial values failed")
          par(mfrow=c(3,1))
          plot(mcmcsamplesD[,"sd_tau"], main= "intial values failed", type='l')
          plot(mcmcsamplesD[,"gamma"], main= "intial values failed", type='l')
          plot(mcmcsamplesD[,"theta"], main= "intial values failed", type='l')
        }
      }
      
      samples[[xx]]  <- mcmcsamples
      samplesS[[xx]]  <- mcmcsamplesS
      samplesD[[xx]]  <- mcmcsamplesD
      
      par(mfrow=c(3,1))
      plot(samplesS[[xx]][ , 'theta'], type='l', main = paste("gamma=",the_gamma_vec[xx]), ylim=c(-5,-2)); 
      abline(cloglog(0.02),0, lty=3, col="red")
      plot(samples[[xx]][ , 'theta'], type='l', main = paste("lambda=", the_lambda_vec[hh]), ylim=c(-5,-2))
      abline(cloglog(0.02),0, lty=3, col="red")
      plot(samplesD[[xx]][ , 'theta'], type='l', main = "", ylim=c(-5,-2)); 
      abline(cloglog(0.02),0, lty=3, col="red")
      
      QD[[xx]] <- c(hdi(samplesD[[xx]][ , 'theta'], credMass=0.9), post_mean=mean(samplesD[[xx]][ , 'theta']))
      QS[[xx]] <- c(hdi(samplesS[[xx]][ , 'theta'], credMass=0.9), post_mean=mean(samplesS[[xx]][ , 'theta']))
      Q[[xx]] <- c(hdi(samples[[xx]][ , 'theta'], credMass=0.9), post_mean=mean(samples[[xx]][ , 'theta']))
      
    }
    
    names(QD)<-paste("gamma =", the_gamma_vec, sep="")
    names(QS)<-paste("gamma =", the_gamma_vec, sep="")
    names(Q)<-paste("gamma =", the_gamma_vec, sep="")
    
    HQ[[hh]]<-Q
    HQS[[hh]]<-QS
    HQD[[hh]]<-QD
    
  }
  names(HQ) <- paste("lambda =", the_lambda_vec, sep="")
  names(HQS) <- paste("lambda =", the_lambda_vec, sep="")
  names(HQD) <- paste("lambda =", the_lambda_vec, sep="")
  
  return(list(HQ=HQ, HQS=HQS, HQD=HQD))}


aa <- sim_once()

######################################################

set.seed(999)

nSim <- 200
#the_lambda_vec <- c(1.5, 0.75, 0.5, 0.33);
the_lambda_vec <- c(1/2, 1/6, 1/10);
the_gamma_vec <- round(c(0, lseq(1, 100, 10)))
true_pop_mu <- 500000

resultsmat <- cbind(expand.grid(the_gamma_vec, the_lambda_vec), 0, 0, 0)
colnames(resultsmat) <- c("gamma", "lambda", "coverage", "width", "improvement")
simmat <- list()
for(hh in 1:length(the_lambda_vec)){ simmat[[hh]] <- list()
for(xx in 1:length(the_gamma_vec)){ simmat[[hh]][[xx]] <- matrix(0, nSim, 3) }}


resultsmatS <- cbind(expand.grid(the_gamma_vec, the_lambda_vec), 0, 0)
colnames(resultsmatS) <- c("gamma", "lambda", "coverage", "width")
simmatS <- list()
for(hh in 1:length(the_lambda_vec)){ simmatS[[hh]] <- list()
for(xx in 1:length(the_gamma_vec)){ simmatS[[hh]][[xx]] <- matrix(0, nSim, 3) }}


resultsmatD <- cbind(expand.grid(the_gamma_vec, the_lambda_vec), 0, 0)
colnames(resultsmatD) <- c("gamma", "lambda", "coverage", "width")
simmatD <- list()
for(hh in 1:length(the_lambda_vec)){ simmatD[[hh]] <- list()
for(xx in 1:length(the_gamma_vec)){ simmatD[[hh]][[xx]] <- matrix(0, nSim, 3) }}




for(isim in 1:nSim){
  print(isim)
  aa <- (sim_once(true_pop_mu = true_pop_mu, 
                  the_gamma_vec = the_gamma_vec, 
                  MCMCiter = 100000, 
                  the_lambda_vec = the_lambda_vec))
  
  for(hh in 1:length(the_lambda_vec)){
    for(xx in 1:length(the_gamma_vec)){	
      simmat[[hh]][[xx]][isim,] <- aa$HQ[[hh]][[xx]];
      simmatS[[hh]][[xx]][isim,] <- aa$HQS[[hh]][[xx]];
      simmatD[[hh]][[xx]][isim,] <- aa$HQD[[hh]][[xx]];
      }}
  
  

  if(isim>3){     
    for(hh in 1:length(the_lambda_vec)){
      for(xx in 1:length(the_gamma_vec)){
        
        resultsmat[resultsmat$gamma==the_gamma_vec[xx] & resultsmat$lambda==the_lambda_vec[hh],"coverage"]<-
          mean(apply(simmat[[hh]][[xx]][1:isim,], 1, function(w) c(w[1]<(cloglog(0.02)) & w[2]>=(cloglog(0.02))) ), na.rm=TRUE)
        resultsmat[resultsmat$gamma==the_gamma_vec[xx] & resultsmat$lambda==the_lambda_vec[hh],"width"]<-
          mean(apply(simmat[[hh]][[xx]][1:isim,], 1, function(w) w[2] - w[1]) ,  na.rm=TRUE)
        
        
        resultsmatS[resultsmatS$gamma==the_gamma_vec[xx] & resultsmatS$lambda==the_lambda_vec[hh],"coverage"]<-
          mean(apply(simmatS[[hh]][[xx]][1:isim,], 1, function(w) c(w[1]<(cloglog(0.02)) & w[2]>=(cloglog(0.02))) ), na.rm=TRUE)
        resultsmatS[resultsmatS$gamma==the_gamma_vec[xx] & resultsmatS$lambda==the_lambda_vec[hh],"width"]<-
          mean(apply(simmatS[[hh]][[xx]][1:isim,], 1, function(w) w[2] - w[1]) ,  na.rm=TRUE)
        
        
        resultsmatD[resultsmatD$gamma==the_gamma_vec[xx] & resultsmatD$lambda==the_lambda_vec[hh],"coverage"]<-
          mean(apply(simmatD[[hh]][[xx]][1:isim,], 1, function(w) c(w[1]<(cloglog(0.02)) & w[2]>=(cloglog(0.02))) ), na.rm=TRUE)
        resultsmatD[resultsmatD$gamma==the_gamma_vec[xx] & resultsmatD$lambda==the_lambda_vec[hh],"width"]<-
          mean(apply(simmatD[[hh]][[xx]][1:isim,], 1, function(w) w[2] - w[1]) ,  na.rm=TRUE)
        
        
        resultsmat[resultsmat$gamma==the_gamma_vec[xx] & resultsmat$lambda==the_lambda_vec[hh],"improvement"]<-
        mean(apply(cbind(simmat[[hh]][[xx]],simmatS[[hh]][[xx]])[1:isim,], 1, function(w) {
          ( 
            ((w[1]<(cloglog(0.02)) & w[2]>=(cloglog(0.02))) > (w[4]<(cloglog(0.02)) & w[5]>=(cloglog(0.02)))) | #cov
              (w[1]<(cloglog(0.02)) & w[2]>=(cloglog(0.02))) & (abs(w[1]-w[2]) < abs(w[4]-w[5]))
            )
          }))
        
      }}
 

    print((cbind(resultsmat,resultsmatS,resultsmatD)))
    if(isim/20==round(isim/20)){
     saveRDS(cbind(resultsmat,resultsmatS,resultsmatD), file=gsub(" ","",paste("simstudypro_",date(), ".rds", sep="")))
      }
  }
}


#saveRDS(cbind(resultsmat,resultsmatS,resultsmatD), file=gsub(" ","",paste("ssfinal2_",date(), ".rds", sep="")))

