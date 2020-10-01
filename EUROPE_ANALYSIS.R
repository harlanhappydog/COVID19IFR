### Bayesian adjustment for preferential testing in estimating infection fatality rates, as motivated by the COVID-19 pandemic
## code: Harlan Campbell
## contact: harlan.campbell@stat.ubc.ca
## details: https://arxiv.org/abs/2005.08459
## github: github.com/harlanhappydog/COVID19IFR/

rm(list = ls(all = TRUE))
ls()

############################
# Determine missing packages and load them:
required_packages <- c("rstan", "ggplot2","xtable", "MCMCvis", "rjags", "coda", "MCMCpack", "BiasedUrn", "rvest", "runjags", "MASS", "HDInterval", "gridGraphics", "gridExtra", "forestplot", "latex2exp", "RCurl","rriskDistributions", "tidyverse")
not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]    
if(length(not_installed)) install.packages(not_installed)                                           
suppressWarnings(lapply(required_packages, require, character.only = TRUE))


############################
# create functions:

	invlog <- function(x){exp(x)/(1+exp(x))}

	logit <- function(x){log(x/(1-x))}

	cloglog <- function(x){log(-log(1-x))}

	icloglog <- function(x){1 - exp(-exp(x))}
	
	
	HDIofMCMC <- function( sampleVec , credMass=0.95 ) {
  # Computes highest density interval from a sample of representative values,
  #   estimated as shortest credible interval.
  # Arguments:
  #   sampleVec
  #     is a vector of representative values from a probability distribution.
  #   credMass
  #     is a scalar between 0 and 1, indicating the mass within the credible
  #     interval that is to be estimated.
  # Value:
  #   HDIlim is a vector containing the limits of the HDI
  sortedPts = sort( sampleVec )
  ciIdxInc = ceiling( credMass * length( sortedPts ) )
  nCIs = length( sortedPts ) - ciIdxInc
  ciWidth = rep( 0 , nCIs )
  for ( i in 1:nCIs ) {
    ciWidth[ i ] = sortedPts[ i + ciIdxInc ] - sortedPts[ i ]
  }
  HDImin = sortedPts[ which.min( ciWidth ) ]
  HDImax = sortedPts[ which.min( ciWidth ) + ciIdxInc ]
  HDIlim = c( HDImin , HDImax )
  return( HDIlim )
}

############################
# stan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Set number of iterations
MCMCiter <- 10000


#####  load datasets
EUROPE_DATA <- getURL("https://raw.githubusercontent.com/harlanhappydog/COVID19IFR/master/EUROPE_DATA.csv")

fullEUR <- read.csv(text = EUROPE_DATA)

SERO_DATA <- getURL("https://raw.githubusercontent.com/harlanhappydog/COVID19IFR/master/SERO_DATA.csv")

phi1_all <- read.csv(text = SERO_DATA)

LOCKDOWN_DATA <- getURL("https://raw.githubusercontent.com/harlanhappydog/COVID19IFR/master/LOCKDOWN_DATA.csv")

lockdown <- read.csv(text = LOCKDOWN_DATA)


##### set seed
set.seed(123)

#######
# Printing data tables for manuscript:

printEUR<-na.omit(fullEUR[,c("Location", "date", "T", "CC", "P", "D",  "aged_70_older", "hospital_beds_per_thousand" , "days_since_first_10infections", "days_till_lockdown", "population_density")])

toprint <- printEUR[,c("Location", "date", "T", "CC", "P", "D")]

rownames(toprint)<-1:dim(toprint)[1]

print(xtable(toprint, digits=c(rep(0,6), rep(2,(dim(toprint)[2]-6)), 0) ),include.rownames=T)

toprint2 <- printEUR[,c("Location",  "aged_70_older", "hospital_beds_per_thousand" , "days_since_first_10infections", "days_till_lockdown", "population_density")]

rownames(toprint2)<-1:dim(toprint2)[1]

print(xtable(toprint2, digits=c(rep(0,3), 2, rep(0,3)) ),include.rownames=T)



#################################################
#################################################
### WRITE OUT MODELS IN STAN
#################################################
#################################################	


#################################################
### MODEL WITH ONLY STERO-STUDIES
#################################################

stan_mod_simple <- "

data {
int<lower=0> kprime;
int<lower=1> K;
int<lower=0> D[K];
int<lower=0> CC[K];
int<lower=0> T[K];
int<lower=0> P[K];
real<lower=0> lambda;
}

parameters {
real<lower=0> gamma;
vector<lower=1, upper=1+gamma>[K] phi;
vector[K] cloglogIFR;
vector[K] cloglogIR;
real theta;
real beta;
real predictIFR;
real predictIR;
real<lower=0> sd_sig;
real<lower=0> sd_tau;
}

transformed parameters {
real<lower=0, upper=1> icloglogtheta;
real<lower=0, upper=1> icloglogbeta;
vector<lower=0, upper=1>[K] IR; 	
vector<lower=0, upper=1>[K] IFR;

icloglogtheta =  1 - exp(-exp(theta));
icloglogbeta =   1 - exp(-exp(beta));

IR  = 1 - exp(-exp(cloglogIR));
IFR = 1 - exp(-exp(cloglogIFR));

}

model {
gamma ~ exponential(lambda);

	phi ~  uniform(1, (1+gamma));	
	cloglogIFR ~ normal(theta, sd_tau);
	cloglogIR ~ normal(beta, sd_sig);

predictIFR ~ normal(theta, sd_tau);
predictIR ~ normal(beta, sd_sig);
	
icloglogtheta ~ uniform(0, 1);
icloglogbeta ~ uniform(0, 1);

sd_sig ~ normal(0, 1);
sd_tau ~ normal(0, 0.1);

for (k in 1:kprime){	
	CC[k] ~ binomial(T[k], IR[k]);	}
if(kprime<K){
for (k in (kprime+1):K){	
	CC[k] ~ binomial(T[k], 1- pow((1-IR[k]), phi[k]) );	}}

	D ~ binomial(P, IR .* IFR);		
}

"

OB_sero <- fullEUR[fullEUR$sero==1,]

datalist_sero <- list(kprime=sum(OB_sero$sero), K=dim(OB_sero)[1],  D=round(OB_sero$D), P=round(OB_sero$P), T=round(OB_sero$T), CC=round(OB_sero$CC), lambda = 0.5)


sero_simple <- stan(model_code = stan_mod_simple,
                 iter = MCMCiter, 
                 warmup = round(MCMCiter*0.1),
                 data = datalist_sero,
                 chains = 4,
                 thin = 50,
                 control = list(max_treedepth = 15, adapt_delta = 0.99)
                 )


## diagnostics
rstan::traceplot(sero_simple, par=c("icloglogtheta","icloglogbeta"))
rstan::traceplot(sero_simple, par=c("predictIFR","predictIR"))

## results
100*icloglog(summary(sero_simple)$summary["predictIFR", c("50%", "2.5%", "97.5%")])
100*icloglog(summary(sero_simple)$summary["predictIR", c("50%", "2.5%", "97.5%")])
100*summary(sero_simple)$summary["icloglogtheta", c("50%", "2.5%", "97.5%")]
100*summary(sero_simple)$summary["icloglogbeta", c("50%", "2.5%", "97.5%")]

main_params <- c("icloglogtheta","icloglogbeta", "theta", "beta", "sd_tau", "sd_sig")

## equal-tailed credible intervals
print(summary(sero_simple)$summary[main_params, c("50%", "2.5%", "97.5%")],3)

## HPD intervals:
summary_md_ci <- function(xx){
c(md=summary(sero_simple)$summary[xx,"50%"], lower=HDIofMCMC(unlist(As.mcmc.list(sero_simple, par=xx)))[1], higher=HDIofMCMC(unlist(As.mcmc.list(sero_simple, par=xx)))[2])}

data.frame(param= main_params, round(t(apply(cbind(main_params),1,function(x) summary_md_ci(x))),3))


OB_sero[,"study_names"]<- factor(paste(1:dim(OB_sero)[1], OB_sero[,"Location"], sep="- "), levels=c(paste(1:dim(OB_sero)[1], OB_sero[,"Location"], sep="- ")))


##########################################
# IFR variables
IFRvars <- cbind(apply(cbind(1:dim(OB_sero)[1]), 1, function(ii) paste("IFR[",ii,"]", sep="")))

study_names <- OB_sero[ ,"study_names"]

meta_IFR <- data.frame(study_names, 100*t(apply(IFRvars, 1, function(xx){summary_md_ci(xx)})))

meta_IFR_plus <- rbind(meta_IFR, data.frame(study_names="Overall", rbind(100*summary_md_ci("icloglogtheta"))))

colnames(meta_IFR_plus) <- c("Study","IFR","lower","upper")
meta_IFR_plus<-droplevels(meta_IFR_plus)

a <- ggplot(meta_IFR_plus, aes(x= Study, y=IFR,ymax=upper,ymin=lower,size=2))+ ylab("IFR (%)")

b <- a+geom_pointrange(cex=c(rep(0.45, datalist_sero$K), 1), col=c(rep("tomato", datalist_sero$kprime), rep("darkolivegreen4", datalist_sero$K-datalist_sero$kprime), "black"))

#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b+coord_flip() + scale_x_discrete(limits = rev(levels(meta_IFR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()


IFR_raw <- phi1_all[order(phi1_all$Location), c("IFR_low",  "IFR_high")]

model_predictIFR <- data.frame(IFR_low = 100*icloglog(summary_md_ci("predictIFR")[2]), IFR_high = 100*icloglog(summary_md_ci("predictIFR")[3]))


model_predictIFR <- model_predictIFR+NA
IFR_extra <- rbind(IFR_raw,model_predictIFR)

IFRplot_sero <- e + geom_pointrange(data=cbind(meta_IFR_plus , IFR_extra), aes(x=Study,ymax=IFR_high,ymin=IFR_low), lwd=.82, pch="", cex=0.3, col=c(rep("lightgrey", dim(IFR_raw)[1]),"lightblue"), position = position_nudge(x = -0.35))+ scale_y_continuous(breaks=seq(0,2,0.25))

IFRplot_sero

##########################################


##########################################
# IR variables
IRvars <- cbind(apply(cbind(1:dim(OB_sero)[1]), 1, function(ii) paste("IR[",ii,"]", sep="")))

study_names <- OB_sero[ ,"study_names"]


meta_IR <- data.frame(study_names, 100*t(apply(IRvars, 1, function(xx){summary_md_ci(xx)})))

meta_IR_plus <- rbind(meta_IR, data.frame(study_names="Overall", rbind(100*summary_md_ci("icloglogbeta"))))

colnames(meta_IR_plus) <- c("Study","IR","lower","upper")
meta_IR_plus <-droplevels(meta_IR_plus)
a <- ggplot(meta_IR_plus, aes(x= Study, y=IR,ymax=upper,ymin=lower,size=2)) + ylab("IR (%)")


b <- a+geom_pointrange(cex=c(rep(0.45, datalist_sero$K), 1), col=c(rep("tomato", datalist_sero$kprime), rep("darkolivegreen4", datalist_sero$K-datalist_sero$kprime), "black"))

#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b+coord_flip() + scale_x_discrete(limits = rev(levels(meta_IR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()


IR_raw <- phi1_all[order(phi1_all$Location), c("IR_low",  "IR_high")]


model_predictIR <- data.frame(IR_low = 100*icloglog(summary_md_ci("predictIR")[2]), IR_high = 100*icloglog(summary_md_ci("predictIR")[3]))


model_predictIR<-model_predictIR+NA
IR_extra<-rbind(IR_raw,model_predictIR)

IRplot_sero <- e + geom_pointrange(data=cbind(meta_IR_plus , IR_extra), aes(x=Study,ymax=IR_high,ymin=IR_low), lwd=.82, pch="", cex=0.3, col=c(rep("lightgrey", dim(IR_raw)[1]),"lightblue"), position = position_nudge(x = -0.35))+labs(x="")

IRplot_sero

##########################################
IFRplot_sero_noyaxis<- IFRplot_sero + theme(axis.title.y = element_blank(),axis.text.y = element_blank())

## Figure 2:
sero_plot <- grid.arrange(IRplot_sero, IFRplot_sero_noyaxis, ncol=2)

# note: requires large plotting area
sero_plot


##########################################

#################################################
### FIT MODEL WITH ALL DATA
#################################################

stan_mod_cov <- "

data {
int<lower=0> kprime;
int<lower=1> K;
int<lower=0> Q;
int<lower=0> H;
int<lower=0> D[K];
int<lower=0> CC[K];
int<lower=0> T[K];
int<lower=0> P[K];
real<lower=0> lambda;
matrix[K,Q] Z;
matrix[K,H] X;
}

parameters {
real<lower=0> gamma;
vector<lower=1, upper=1+gamma>[K] phi;
vector[K] cloglogIFR;
vector[K] cloglogIR;
real theta;
real beta;
vector[Q] thetavec;
vector[H] betavec;

real predictIFR;
real predictIR;
real<lower=0> sd_sig;
real<lower=0> sd_tau;
}

transformed parameters {
real<lower=0, upper=1> icloglogtheta;
real<lower=0, upper=1> icloglogbeta;
vector<lower=0, upper=1>[K] IR; 	
vector<lower=0, upper=1>[K] IFR;

icloglogtheta =  1 - exp(-exp(theta));
icloglogbeta =   1 - exp(-exp(beta));

IR = 1 - exp(-exp(cloglogIR));
IFR = 1 - exp(-exp(cloglogIFR));
}

model {
thetavec ~ normal(0, 1);
betavec ~ normal(0, 1);

gamma ~ exponential(lambda);
phi ~  uniform(1, (1+gamma));	

cloglogIFR ~ normal(Z*thetavec + theta, sd_tau);
cloglogIR ~ normal(X*betavec + beta, sd_sig);

icloglogtheta ~ uniform(0, 1);
icloglogbeta ~ uniform(0, 1);

predictIFR ~ normal(theta, sd_tau);
predictIR ~ normal(beta, sd_sig);
	
sd_sig ~ normal(0, 1);
sd_tau ~ normal(0, 0.1);

for (k in 1:kprime){	
	CC[k] ~ binomial(T[k], IR[k]);	}
if(kprime<K){
for (k in (kprime+1):K){	
	CC[k] ~ binomial(T[k], 1- pow((1-IR[k]), phi[k]) );	}}
	
D ~ binomial(P, IR .* IFR);	
	
}

"

OB <- fullEUR[!rowSums(is.na(fullEUR[,c("scale_days_till_lockdown", "scale_population_density", "scale_days_since_first_10infections", "scale_aged_70_older", "scale_hospital_beds_per_thousand")])),]

datalist <- list(kprime = sum(OB$sero), K = dim(OB)[1], H = 3, Q = 2, D = round(OB$D), P = round(OB$P), T = round(OB$T), CC = round(OB$CC), lambda = 0.5, X = cbind(c(OB$scale_days_since_first_10infections), c(OB$scale_days_till_lockdown), c(OB$scale_population_density)),  Z = cbind(c(OB$scale_aged_70_older), c(OB$scale_hospital_beds_per_thousand)))


MAAD <- stan(model_code = stan_mod_cov,
                 iter = MCMCiter, 
                 warmup = round(MCMCiter*0.1),
                 data = datalist,
                 chains = 4,
                 thin = 50,
                 control = list(max_treedepth = 15, adapt_delta = 0.99)
                 )


## diagnostics
rstan::traceplot(MAAD, par=c("icloglogtheta","icloglogbeta"))
rstan::traceplot(MAAD, par=c("thetavec[1]","thetavec[2]"))
rstan::traceplot(MAAD, par=c("betavec[1]","betavec[2]"))
rstan::traceplot(MAAD, par=c("predictIR","predictIFR"))
rstan::traceplot(MAAD, par=c("gamma"))

## results
100*summary(MAAD)$summary["icloglogtheta", c("50%", "2.5%", "97.5%")]
100*summary(MAAD)$summary["icloglogbeta", c("50%", "2.5%", "97.5%")]


main_params <- c("icloglogtheta","icloglogbeta", "theta", "beta", "thetavec[1]", "thetavec[2]",  "betavec[1]", "betavec[2]","betavec[3]", "sd_tau", "sd_sig", "gamma")

## equal-tailed credible intervals
print(summary(MAAD)$summary[main_params, c("50%", "2.5%", "97.5%")],3)


## HPD intervals:
summary_md_ci_MAAD <- function(xx){
c(md=summary(MAAD)$summary[xx,"50%"], lower=HDIofMCMC(unlist(As.mcmc.list(MAAD, par=xx)))[1], higher=HDIofMCMC(unlist(As.mcmc.list(MAAD, par=xx)))[2])}

data.frame(param=c(main_params),round(t(apply(cbind(c(main_params)),1,function(x) summary_md_ci_MAAD(x))),3))


OB[,"study_names"]<- factor(paste(1:dim(OB)[1], OB[,"Location"], sep="- "), levels=c(paste(1:dim(OB)[1], OB[,"Location"], sep="- ")))



##########################################
# IFR variables
IFRvars <- cbind(apply(cbind(1:dim(OB)[1]), 1, function(ii) paste("IFR[",ii,"]", sep="")))

study_names <- OB[ ,"study_names"]

meta_IFR <- data.frame(study_names, 100*t(apply(IFRvars, 1, function(xx){summary_md_ci_MAAD(xx)})))


meta_IFR_plus <- rbind(meta_IFR, data.frame(study_names="Overall", rbind(100*summary_md_ci_MAAD("icloglogtheta"))))

colnames(meta_IFR_plus) <- c("Study", "IFR","lower","upper")

a <- ggplot(meta_IFR_plus, aes(x= Study, y=IFR,ymax=upper,ymin=lower,size=2))+ ylab("IFR (%)")


b <- a+geom_pointrange(cex=c(rep(0.45, datalist$K), 0.5), col=c(rep("tomato", datalist$kprime), rep("darkolivegreen4", datalist$K-datalist$kprime), "black"))

#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b+coord_flip() + scale_x_discrete(limits = rev(levels(meta_IFR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()


IFR_raw <- phi1_all[order(phi1_all$Location), c("IFR_low",  "IFR_high")]

NAmat<-matrix(rep(c(NA,NA), dim(meta_IFR)[1]-dim(IFR_raw)[1]),,2)
colnames(NAmat)<-c("IFR_low",  "IFR_high")

model_predictIFR <- data.frame(IFR_low = 100*icloglog(summary_md_ci_MAAD("predictIFR")[2]), IFR_high = 100*icloglog(summary_md_ci_MAAD("predictIFR")[3]))


model_predictIFR <-model_predictIFR+NA
IFR_extra<-rbind(IFR_raw, NAmat, model_predictIFR)


IFRplot<-e+ geom_pointrange(data=cbind(meta_IFR_plus , IFR_extra), aes(x=Study,ymax=IFR_high,ymin=IFR_low), lwd=.82, pch="", cex=0.3, col=c(rep("lightgrey", dim(rbind(IFR_raw,NAmat))[1]),"lightblue"), position = position_nudge(x = -0.35))
IFRplot


##########################################


##########################################
# IR variables
IRvars <- cbind(apply(cbind(1:dim(OB)[1]), 1, function(ii) paste("IR[",ii,"]", sep="")))

study_names <- OB[,"study_names"]

meta_IR <- data.frame(study_names, 100*t(apply(IRvars, 1, function(xx){summary_md_ci_MAAD(xx)})))

meta_IR_plus <- rbind(meta_IR, data.frame(study_names="Overall", rbind(100*summary_md_ci_MAAD("icloglogbeta"))))

colnames(meta_IR_plus) <- c("Study","IR","lower","upper")
meta_IR_plus <-droplevels(meta_IR_plus)


a <- ggplot(meta_IR_plus, aes(x= Study, y=IR,ymax=upper,ymin=lower,size=2)) + ylab("IR (%)")

b <- a+geom_pointrange(cex=c(rep(0.45, datalist$K), 0.5), col=c(rep("tomato", datalist$kprime), rep("darkolivegreen4", datalist$K-datalist$kprime), "black"))

#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b+coord_flip() + scale_x_discrete(limits = rev(levels(meta_IR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()

IRplot <- e 
IRplot

IR_raw <- phi1_all[order(phi1_all$Location), c("IR_low",  "IR_high")]

NAmat<-matrix(rep(c(NA,NA),dim(meta_IR_plus)[1]-dim(IR_raw)[1]),,2)
colnames(NAmat)<-c("IR_low",  "IR_high")
IRplot <- e + geom_pointrange(data=cbind(meta_IR_plus , rbind(IR_raw, NAmat)), aes(x=Study,ymax=IR_high,ymin=IR_low), lwd=.82, pch="", cex=0.3, col="lightgrey", position = position_nudge(x = -0.35))

IRplot



IR_raw <- phi1_all[order(phi1_all$Location), c("IR_low", "IR_high")]

NAmat<-matrix(rep(c(NA,NA), dim(meta_IR)[1]-dim(IR_raw)[1]),,2)
colnames(NAmat)<-c("IR_low",  "IR_high")

model_predictIR <- data.frame(IR_low = 100*icloglog(summary_md_ci_MAAD("predictIR")[2]), IR_high = 100*icloglog(summary_md_ci_MAAD("predictIR")[3]))

model_predictIR <- model_predictIR+NA
IR_extra<-rbind(IR_raw, NAmat, model_predictIR)


IRplot<-e+ geom_pointrange(data=cbind(meta_IR_plus , IR_extra), aes(x=Study,ymax=IR_high,ymin=IR_low), lwd=.82, pch="", cex=0.3, col=c(rep("lightgrey", dim(rbind(IR_raw,NAmat))[1]),"lightblue"), position = position_nudge(x = -0.35))+labs(x="")
IRplot



##########################################
IFRplot_noyaxis<- IFRplot + theme(axis.title.y = element_blank(),axis.text.y = element_blank())

# Figure 3
grid.arrange(IRplot, IFRplot_noyaxis, ncol=2)


diff(100*HDIofMCMC(unlist(As.mcmc.list(sero_simple, par="predictIFR"))))
diff(100*HDIofMCMC(unlist(As.mcmc.list(MAAD, par="predictIFR"))))
            
diff(100*HDIofMCMC(unlist(As.mcmc.list(sero_simple, par="icloglogtheta"))))
diff(100*HDIofMCMC(unlist(As.mcmc.list(MAAD, par="icloglogtheta"))))


round(100*HDIofMCMC(unlist(As.mcmc.list(sero_simple, par="icloglogtheta"))),3)
round(100*HDIofMCMC(unlist(As.mcmc.list(MAAD, par="icloglogtheta"))),3)
            
      
round(100*summary(sero_simple)$summary["icloglogtheta", c("50%")],3)
round(100*HDIofMCMC(unlist(As.mcmc.list(sero_simple, par="icloglogtheta"))),3)

round(100*summary(MAAD)$summary["icloglogtheta", c("50%")],3)
round(100*HDIofMCMC(unlist(As.mcmc.list(MAAD, par="icloglogtheta"))),3)


       
##########################################
# Scaterplot of phi_k vs. H2 index
##########################################

#ww<-getURL("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv")
#lockdown<-read.csv(text=ww)
#write.csv(lockdown, '~/Desktop/UBC/RECODID_ZIKV/COVID/jasa_template/code/LOCKDOWN_DATA.csv')
kprime<-5

testingdays <- lockdown [as.Date(as.character(lockdown[,"Date"]),  "%Y%m%d")> "2020-02-01" & as.Date(as.character(lockdown[,"Date"]),  "%Y%m%d")< "2020-04-01" & lockdown[,"CountryCode"]%in%fullEUR$iso_code,c("Date", "CountryCode","H2_Testing.policy"),]

testingdays <-unique(testingdays)

library("tidyverse")
testingdays_wide <- testingdays %>% spread("Date","H2_Testing.policy")

testing_avg <- data.frame(CountryCode=testingdays_wide[,1],testing_avg=rowMeans(testingdays_wide[,-1]))


phivars<-cbind(apply(cbind(1:dim(OB)[1]), 1, function(ii) paste("phi[",ii,"]", sep="")))
phivars<-phivars[-c(1:kprime)]

phi_est<-data.frame(CountryCode =OB$iso_code[-c(1:kprime)], phi= summary(MAAD)$summary[phivars,"mean"])

compare_testing<-merge(testing_avg,phi_est, by="CountryCode", all=TRUE)

compare_testing2<-merge(compare_testing, fullEUR[fullEUR[,"sero"]==0,], by.x="CountryCode", by.y= "iso_code")

gg <- ggplot(compare_testing2,
       mapping= aes(x= phi, y = testing_avg))+
  geom_point(aes(size=P), shape = 21, color="white", fill = "lightblue")+

  ggrepel::geom_text_repel(aes(label = Location),
                           size = 3, segment.color = NA,
                           point.padding = unit(0.1, "lines")) +
  theme_classic() +

  # This scales area to size (not radius), specifies max size, and hides legend
  scale_size_area(max_size = 20, guide = FALSE)



legend_bubbles <- data.frame(
  label = c("5m", "15m", "40m"),
  size  = c(5E6, 15E6, 40E6)
) %>%
mutate(radius = sqrt(size / pi))   


gg + scale_y_continuous(limits = c(0, 1.75))+geom_point(data = legend_bubbles,
             #  The "radius/301" was trial and error. Better way?
             aes(x = 1.65, y = 0.0125 + (radius/345)/140, size = size),
             shape = 21, color = "black", fill = NA) +
  geom_text(data = legend_bubbles, size = 2.5,
            aes(x = 1.65, y = 0.027 + (2.15 * radius/351)/140, label = label)) +
  annotate("text", x = 1.65, y = 0.24, label = "Population", fontface = "bold") + labs(y= "H2 index", x=expression(phi[k]), size=3) + theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) 
        #+geom_smooth(span=2, se=FALSE, mapping = aes(weight = P))
