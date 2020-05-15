############################################################
#### Bayesian adjustment for preferential testing ####
#### in estimating the COVID-19 infection fatality rate ####
#### code written by Harlan Campbell ####
#### contact harlan.campbell@stat.ubc.ca ####

#### Code to reproduce illustrative example:

# Determine missing packages and load them:
required_packages <- c("MCMCvis", "rjags", "coda", "MCMCpack", "BiasedUrn", "rvest", "runjags", "MASS", "HDInterval", "gridGraphics", "gridExtra", "forestplot", "latex2exp")
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

sd(OD$df$IFR)

#install.packages("quadprog")
library(quadprog)

minvar <- function(mn, theta.lo, theta.hi) {

  ### first check constraints consistent, return NA if not
  ans <- NA
  if ( (mean(theta.lo) < mn) && (mn < mean(theta.hi)) ) {
  
    ### reparam as theta[i] = theta.lo[j] + x[i]*(theta.hi[j]-theta.lo[j])
  
    J <- length(theta.lo)
    wdth <- theta.hi - theta.lo
  
    ### linear equality/inequality constraints Amat%*%x >= bvec
    Amat <- matrix(NA, J, 2*J+1); bvec <- rep(NA, 2*J+1)
  
    ### one equality constraint
    Amat[,1] <- wdth/J;  bvec[1] <- mn - mean(theta.lo)
  
    ### 2J inequality constraints
    Amat[,2:(J+1)] <- diag(J); bvec[2:(J+1)] <- rep(0,J)
    Amat[,(J+2):(2*J+1)] <- -diag(J); bvec[(J+2):(2*J+1)] <- rep(-1,J)

    ### minimization target is 0.5*sum((th-mn)^2)  
    ### expressed as .5 t(x)%*%Dmat%*%x - t(dvec)%*%x
    Dmat <- diag(wdth^2)
    dvec <- -wdth*theta.lo
  
    ans <- solve.QP(Dmat, dvec, Amat, bvec, meq=1)
  }
ans
}

IDint <- function(bnd, a, b, phi.lo, phi.hi, theta.tr, tol=0.001) {

### bnd is the upper bound on SD(IFR)
  
### a,b, phi.lo, phi.hi all length I vectors,
  
### theta.tr (or IFR.overbar.true) useful input,
### since need to start the numerical search at a value
### known to be in the ID interval  
  
  theta.lo <- a/(1-(1-b)^(1/phi.lo))
  theta.hi <- a/(1-(1-b)^(1/phi.hi))  
  
  upr <- theta.tr; flg <- T
  while (flg) {
    upr <- upr + tol
    tmp <- minvar(upr, theta.lo, theta.hi)
    flg <- F
    if (!is.na(tmp[[1]][1])) {
      tmp2 <- theta.lo + tmp$sol*(theta.hi-theta.lo)
      if (sqrt(var(tmp2)) < bnd ) {
        flg <- T
      } 
    }
  }  
  upr <- upr - tol

  lwr <- theta.tr; flg <- T
  while (flg) {
    lwr <- lwr - tol
    tmp <- minvar(lwr, theta.lo, theta.hi)
    flg <- F
    if (!is.na(tmp[[1]][1])) {
      tmp2 <- theta.lo + tmp$sol*(theta.hi-theta.lo)
      if (sqrt(var(tmp2)) < bnd ) {
        flg <- T
      }
    }
  }  
  lwr <- lwr + tol

c(lwr,upr)  
}




### number of jurisdictions
k <- 12

### average IFR
theta <- 0.02
tau_max <- 0.002
sd(OD$df$IFR)


### various IRs
lambda <- c(140, 206, 303, 190, 159, 132, 137, 526, 212, 160, 164, 245)/1000
phi_list<-list()
phi_list[[3]] <- seq(1,23, length.out=12)
phi_list[[2]] <- seq(1,12)
phi_list[[1]] <- seq(1,5, length.out=12)
gamma_vec <- c(4,11,22)

layout(matrix(c(1,2,3,4,5,6),3,2))
par(las=c(1))
for(ii in 1:3){

  ### mix of preferentialities
  phi<-phi_list[[ii]]
  a <- theta*lambda
  b <- 1-(1-lambda)^phi
  
  phi.lo <- rep(1,k)
  phi.hi <- rep(40,k)
  
  
  if(ii==1){
    plot(-1,-1, xlim=c(1,k),ylim=c(0,.25), xlab="k",     ylab="IFR", main= expression(paste(gamma," = 4,   ", bar(tau), " = 0")),xaxt ="n");axis(1, at=1:12)}
  if(ii==2){
    plot(-1,-1, xlim=c(1,k),ylim=c(0,.25), xlab="k",     ylab="IFR", main= expression(paste(gamma," = 11,   ", bar(tau), " = 0")),xaxt ="n");axis(1, at=1:12)}
  if(ii==3){
    plot(-1,-1, xlim=c(1,k),ylim=c(0,.25), xlab="k",     ylab="IFR", main= expression(paste(gamma," = 22,   ", bar(tau), " = 0")),xaxt ="n");axis(1, at=1:12)}
  
  theta.lo <- a/(1-(1-b)^(1/phi.lo))
  theta.hi <- a/(1-(1-b)^(1/phi.hi))
  
  rect(1,max(theta.lo),k,min(theta.hi), col="lightgreen", border=NA)
  print(ii)
  print(round(c(max(theta.lo),min(theta.hi)),4))
  for (i in 1:k) {
    points(rep(i,2),c(theta.lo[i], theta.hi[i]), 
           type="l", lwd=2, lend='butt')       
  }
 
}


for(ii in 1:3){
  
  ### mix of preferentialities
  phi<-phi_list[[ii]]
  a <- theta*lambda
  b <- 1-(1-lambda)^phi
  
  phi.lo <- rep(1,k)
  phi.hi <- rep(40,k)
  

  theta.lo <- a/(1-(1-b)^(1/phi.lo))
  theta.hi <- a/(1-(1-b)^(1/phi.hi))
  
  
  if(ii==1){
    plot(-1,-1, xlim=c(1,k),ylim=c(0,.25), xlab="k",     ylab="IFR", main= expression(paste(gamma," = 4,   ", bar(tau), " = 0.002")),xaxt ="n");axis(1, at=1:12)}
  if(ii==2){
    plot(-1,-1, xlim=c(1,k),ylim=c(0,.25), xlab="k",     ylab="IFR", main= expression(paste(gamma," = 11,   ", bar(tau), " = 0.002")),xaxt ="n");axis(1, at=1:12)}
  if(ii==3){
    plot(-1,-1, xlim=c(1,k),ylim=c(0,.25), xlab="k",     ylab="IFR", main= expression(paste(gamma," = 22,   ", bar(tau), " = 0.002")),xaxt ="n");axis(1, at=1:12)}
  
  gblim<-IDint(tau_max, a=theta*lambda, b= 1-(1-lambda)^phi, phi.lo=rep(1,k), phi.hi=rep(40,k), theta.tr=0.02, tol=0.0001)  
  
  rect(1,gblim[1],k,gblim[2], col="lightgreen", border=NA)
  print(paste(ii," tau_bar=",tau_max))
  print(round(c(gblim[1],gblim[2]),4))
  
  for (i in 1:k) {
    points(rep(i,2),c(theta.lo[i], theta.hi[i]), 
           type="l", lwd=2, lend='butt')       
  }
}  




