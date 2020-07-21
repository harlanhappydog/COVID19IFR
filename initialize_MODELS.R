############
# NIMBLE MODEL "D"
############

modelCodeD <- nimbleCode( {
  
  #priors
  for (k in (kprime+1):K){  	phi[k] ~ dunif(1, 1+gamma);}
  
  gamma  ~ dexp(lambda);
  icloglog_theta ~ dunif(0, 1); 
  icloglog_beta0 ~ dunif(0, 1);
  theta <- log(-log(1-icloglog_theta));
  beta0 <- log(-log(1-icloglog_beta0));   
  sd_sig   ~ dnorm(0, sd=1) ;
  sd_tau   ~ dnorm(0, sd=0.1) ;
  constraint_data1 ~ dconstraint( sd_sig > 0 )
  constraint_data2 ~ dconstraint( sd_tau > 0 )
  
  #likelihood 
  for (k in (kprime+1):K){
    cloglog_IFR[k]              ~ dnorm(theta, sd=sd_tau)
    cloglog_infectionrate[k]    ~ dnorm(beta0, sd=sd_sig)
    cloglog(infectionrate[k])   <- cloglog_infectionrate[k]
    cloglog(IFR[k])             <- cloglog_IFR[k]
    
    confirmed_cases[k] ~ dbin(1-(1-infectionrate[k])^(phi[k]), tests[k]);
    deaths[k] ~ dbin(IFR[k] * infectionrate[k], population[k])
  }
})

############

OB_data <- simulate_data(K = K_sim,
                         kprime = kprime_sim,
                         theta = cloglog(0.02), 
                         tausquared = 0.005, 
                         beta0 = cloglog(0.2), 
                         sigmasquared = 0.25, 
                         pop_mu = 10000, 
                         testingrate_lower = 0.001, 
                         testingrate_upper = 0.01, 
                         gamma_vec = c(0,5,20,50))


# Assembliong all data model:
lambda_now <- 0.5
consts <- list(K = length(OB_data$df$kvec), kprime = kprime_sim)
data <- list(confirmed_cases = OB_data$CC[[1]], 
             deaths = OB_data$df$deaths, 
             population =  OB_data$df$population, 
             tests = OB_data$df$tests, 
             lambda = lambda_now)

# Reasonable initial values:
cloglog_infectionrate_init <- runif(consts$K, 0, 0.5)
cloglog_IFR_init <- runif(consts$K, 0, 0.1)
gamma_intial <- rexp(1, lambda_now)
inits <- list(cloglog_IFR = cloglog_IFR_init,
              cloglog_infectionrate = cloglog_infectionrate_init,
              gamma = gamma_intial,
              phi = c(rep(-99,consts$kprime), runif(consts$K-consts$kprime, 1, 1+gamma_intial)),
              icloglog_theta = mean(cloglog_IFR_init),
              icloglog_beta0 = mean(cloglog_infectionrate_init),
              sd_tau = abs(rnorm(1, sd=0.1)),
              sd_sig = abs(rnorm(1, sd=1)))

modelD <- nimbleModel(code = modelCodeD,
                     name = 'covidD',
                     constants = consts,
                     data = data,
                     inits = inits)

cloglog_IFR_nodes <- modelD$expandNodeNames('cloglog_IFR')
cloglog_infectionrate_nodes <- modelD$expandNodeNames('cloglog_infectionrate')
phi_nodes <- modelD$expandNodeNames('phi')

####

MCMCconf2D <- configureMCMC(modelD, monitors = c('phi','gamma','theta','beta0','sd_tau','sd_sig','cloglog_IFR','cloglog_infectionrate'))
source('~/Desktop/UBC/RECODID_ZIKV/COVID/Rcode/samplers.R', chdir = TRUE)
for(i in seq_along(phi_nodes)) {
  ## We could remove some default samplers, but I'll leave them in for now.
  MCMCconf2D$addSampler(type = 'sampler_infrate_phi',
                       target = c(cloglog_infectionrate_nodes[i],
                                  phi_nodes[i]),
                       silent = TRUE)
  MCMCconf2D$addSampler(type = 'sampler_infrate_IFR',
                       target = c(cloglog_infectionrate_nodes[i],
                                  cloglog_IFR_nodes[i]),
                       silent = TRUE)
}
#MCMCconf2$printSamplers()

MCMC2D    <- buildMCMC(MCMCconf2D)
CmodelD   <- compileNimble(modelD)
CMCMCsD   <- compileNimble(MCMC2D, project = modelD)
CMCMC2D   <- CMCMCsD
samples2D <- suppressWarnings(runMCMC(CMCMC2D, thin = 10, nburnin = 100, niter = 10000))
plot(samples2D[,"beta0"], type='l')
plot(samples2D[,"theta"], type='l')




############
# NIMBLE MODEL
############

modelCode <- nimbleCode( {
  
  #priors
  for (k in 1:kprime){ 		phi[k] ~ dunif(1,1.0000001);}
  for (k in (kprime+1):K){  	phi[k] ~ dunif(1, gamma+1);}
  
  gamma  ~ dexp(lambda);
  icloglog_theta ~ dunif(0, 1); 
  icloglog_beta0 ~ dunif(0, 1);
  theta <- log(-log(1-icloglog_theta));
  beta0 <- log(-log(1-icloglog_beta0));   
  sd_sig   ~ dnorm(0, sd=1) ;
  sd_tau   ~ dnorm(0, sd=0.1) ;
  constraint_data1 ~ dconstraint( sd_sig > 0 )
  constraint_data2 ~ dconstraint( sd_tau > 0 )
  
  #likelihood 
  for (k in 1:K){
    cloglog_IFR[k]              ~ dnorm(theta, sd=sd_tau)
    cloglog_infectionrate[k]    ~ dnorm(beta0, sd=sd_sig)
    cloglog(infectionrate[k])   <- cloglog_infectionrate[k]
    cloglog(IFR[k])             <- cloglog_IFR[k]
    
    confirmed_cases[k] ~ dbin(1-(1-infectionrate[k])^(phi[k]), tests[k]);
    deaths[k] ~ dbin(IFR[k] * infectionrate[k], population[k])
  }
})

############




# Assembliong all data model:
lambda_now <- 0.5
consts <- list(K = length(OB_data$df$kvec), kprime = kprime_sim)
data <- list(confirmed_cases = OB_data$CC[[1]], 
             deaths = OB_data$df$deaths, 
             population =  OB_data$df$population, 
             tests = OB_data$df$tests, 
             lambda = lambda_now)

# Reasonable initial values:
cloglog_infectionrate_init <- runif(consts$K, 0, 0.5)
cloglog_IFR_init <- runif(consts$K, 0, 0.1)
gamma_intial <- rexp(1, lambda_now)
inits <- list(cloglog_IFR = cloglog_IFR_init,
              cloglog_infectionrate = cloglog_infectionrate_init,
              gamma = gamma_intial,
              phi = c(rep(1,consts$kprime), runif(consts$K-consts$kprime, 1, 1+gamma_intial)),
              icloglog_theta = mean(cloglog_IFR_init),
              icloglog_beta0 = mean(cloglog_infectionrate_init),
              sd_tau = abs(rnorm(1, sd=0.1)),
              sd_sig = abs(rnorm(1, sd=1)))


model <- nimbleModel(code = modelCode,
                     name = 'covid',
                     constants = consts,
                     data = data,
                     inits = inits)

cloglog_IFR_nodes <- model$expandNodeNames('cloglog_IFR')
cloglog_infectionrate_nodes <- model$expandNodeNames('cloglog_infectionrate')
phi_nodes <- model$expandNodeNames('phi')

####

MCMCconf2 <- configureMCMC(model, monitors = c('phi','gamma','theta','beta0','sd_tau','sd_sig','cloglog_IFR','cloglog_infectionrate'))
source('~/Desktop/UBC/RECODID_ZIKV/COVID/Rcode/samplers.R', chdir = TRUE)
for(i in seq_along(phi_nodes)) {
  ## We could remove some default samplers, but I'll leave them in for now.
  MCMCconf2$addSampler(type = 'sampler_infrate_phi',
                       target = c(cloglog_infectionrate_nodes[i],
                                  phi_nodes[i]),
                       silent = TRUE)
  MCMCconf2$addSampler(type = 'sampler_infrate_IFR',
                       target = c(cloglog_infectionrate_nodes[i],
                                  cloglog_IFR_nodes[i]),
                       silent = TRUE)
}
#MCMCconf2$printSamplers()

MCMC2    <- buildMCMC(MCMCconf2)
Cmodel   <- compileNimble(model)
CMCMCs   <- compileNimble(MCMC2, project = model)
CMCMC2   <- CMCMCs
samples2 <- suppressWarnings(runMCMC(CMCMC2, thin = 10, nburnin = 100, niter = 10000))
plot(samples2[,"beta0"], type='l')
plot(samples2[,"theta"], type='l')



############
# NIMBLE MODEL "S"
############

modelCodeS <- nimbleCode( {
  
  #priors
  for (k in 1:kprime){ 		phi[k] ~ dunif(1,1.0000001);}
  
  icloglog_theta ~ dunif(0, 1); 
  icloglog_beta0 ~ dunif(0, 1);
  theta <- log(-log(1-icloglog_theta));
  beta0 <- log(-log(1-icloglog_beta0));   
  sd_sig   ~ dnorm(0, sd=1) ;
  sd_tau   ~ dnorm(0, sd=0.1) ;
  constraint_data1 ~ dconstraint( sd_sig > 0 )
  constraint_data2 ~ dconstraint( sd_tau > 0 )
  
    
  #likelihood 
  for (k in 1:kprime){
    cloglog_IFR[k]              ~ dnorm(theta, sd=sd_tau)
    cloglog_infectionrate[k]    ~ dnorm(beta0, sd=sd_sig)
    cloglog(infectionrate[k])   <- cloglog_infectionrate[k]
    cloglog(IFR[k])             <- cloglog_IFR[k]
    
    confirmed_cases[k] ~ dbin(1-(1-infectionrate[k])^(phi[k]), tests[k]);
    deaths[k] ~ dbin(IFR[k] * infectionrate[k], population[k])
  }
})

############


# Assembliong all data model:
consts <- list(K = length(OB_data$df$kvec), kprime = kprime_sim)
data <- list(confirmed_cases = OB_data$CC[[1]], 
             deaths = OB_data$df$deaths, 
             population =  OB_data$df$population, 
             tests = OB_data$df$tests)

# Reasonable initial values:
cloglog_infectionrate_init <- runif(consts$K, 0, 0.5)
cloglog_IFR_init <- runif(consts$K, 0, 0.1)
inits <- list(cloglog_IFR = cloglog_IFR_init,
              cloglog_infectionrate = cloglog_infectionrate_init,
              icloglog_theta = mean(cloglog_IFR_init),
              icloglog_beta0 = mean(cloglog_infectionrate_init),
              sd_tau = abs(rnorm(1, sd=0.1)),
              sd_sig = abs(rnorm(1, sd=1)))


modelS <- nimbleModel(code = modelCodeS,
                      name = 'covidS',
                      constants = consts,
                      data = data,
                      inits = inits)

cloglog_IFR_nodes <- modelS$expandNodeNames('cloglog_IFR')
cloglog_infectionrate_nodes <- modelS$expandNodeNames('cloglog_infectionrate')
phi_nodes <- modelS$expandNodeNames('phi')

####

MCMCconf2S <- configureMCMC(modelS, monitors = c('phi','theta','beta0','sd_tau','sd_sig','cloglog_IFR','cloglog_infectionrate'))
source('~/Desktop/UBC/RECODID_ZIKV/COVID/Rcode/samplers.R', chdir = TRUE)
for(i in seq_along(phi_nodes)) {
  ## We could remove some default samplers, but I'll leave them in for now.
  MCMCconf2S$addSampler(type = 'sampler_infrate_phi',
                        target = c(cloglog_infectionrate_nodes[i],
                                   phi_nodes[i]),
                        silent = TRUE)
  MCMCconf2S$addSampler(type = 'sampler_infrate_IFR',
                        target = c(cloglog_infectionrate_nodes[i],
                                   cloglog_IFR_nodes[i]),
                        silent = TRUE)
}
#MCMCconf2S$printSamplers()

MCMC2S    <- buildMCMC(MCMCconf2S)
CmodelS   <- compileNimble(modelS)
CMCMCsS   <- compileNimble(MCMC2S, project = modelS)
CMCMC2S   <- CMCMCsS
samples2S <- suppressWarnings(runMCMC(CMCMC2S, thin = 10, nburnin = 100, niter = 10000))
plot(samples2S[,"beta0"], type='l')
plot(samples2S[,"theta"], type='l')

