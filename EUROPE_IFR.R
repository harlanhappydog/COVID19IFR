
############################
# Determine missing packages and load them:
required_packages <- c("rstan", "ggplot2","xtable", "MCMCvis", "rjags", "coda", "MCMCpack", "BiasedUrn", "rvest", "runjags", "MASS", "HDInterval", "gridGraphics", "gridExtra", "forestplot", "latex2exp","RCurl","rriskDistributions")
not_installed <- required_packages[!(required_packages %in% installed.packages()[ , "Package"])]    
if(length(not_installed)) install.packages(not_installed)                                           
suppressWarnings(lapply(required_packages, require, character.only = TRUE))


############################
# create functions:
	invlog <- function(x){exp(x)/(1+exp(x))}
	logit <- function(x){log(x/(1-x))}
	cloglog <- function(x){log(-log(1-x))}
	icloglog <- function(x){1 - exp(-exp(x))}
	
	
HDIofMCMC = function( sampleVec , credMass=0.95 ) {
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


#################################################
#################################################
### BUILD THE DATASET
#################################################
#################################################	

x <- getURL("https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv")
owid <- read.csv(text = x)
#write.csv(owid, file= '~/Desktop/UBC/RECODID_ZIKV/COVID/Rcode/owid_july16_2020.csv')


# only Europe, and only first week of May:
EUR_may1_all <- (owid[owid[,"continent"]%in%"Europe" & owid[,"date"]%in%c("2020-05-01","2020-05-02","2020-05-03","2020-05-04","2020-05-05","2020-05-06","2020-05-07") ,])


# only need certain variables:
EUR_may1 <- EUR_may1_all[,c("iso_code", "location", "date", "total_cases", "total_tests", "total_cases_per_million", "aged_70_older", "hospital_beds_per_thousand")]

# only need earliest data available during first week of May:
country <- unique(EUR_may1_all[,"iso_code"])
EUR_j <- NULL
for(j in 1:length(unique(EUR_may1_all[,"iso_code"]))){
EUR_j <- rbind(EUR_j, na.omit(EUR_may1[EUR_may1[,"iso_code"]%in%country[j],])[1,])
}

EUR_may1 <- na.omit(EUR_j)
dim(EUR_may1)

# Note that France, Sweden and others are missing due to lack of testing data:
EUR_may1_all[EUR_may1_all[,"location"]%in%c("France", "Sweden"),]


# Add population data:
EUR_may1[,"Population"] <- (EUR_may1[,"total_cases"]/EUR_may1[,"total_cases_per_million"])*(1000000)

# obtain deaths data from +14 days
EUR_may1[,"total_deaths"] <- NULL
country <- unique(EUR_may1_all[,"iso_code"])
for(icountry in country){
	
test_date  <-	EUR_may1[EUR_may1[, "iso_code"]%in%icountry, "date"] 
death_date <- as.Date(test_date) +14

EUR_may1[EUR_may1[, "iso_code"]%in%icountry, "total_deaths"] <- owid[owid[,"continent"]=="Europe"  & owid[, "iso_code"]%in%c(icountry) & owid[,"date"]%in%paste(death_date), c("total_deaths")]
}

EUR_dat<-EUR_may1

# restrict countries to EEA country club + GB + CHE
EEAclub <- c("Austria", "Belgium", "Bulgaria", "Croatia", "Cyprus", "Czech Republic", "Denmark", "Estonia", "Finland", "France", "Germany", "Greece", "Hungary", "Ireland", "Italy", "Latvia", "Lithuania", "Luxembourg", "Malta", "Netherlands", "Poland", "Portugal", "Romania", "Slovakia", "Slovenia", "Spain", "Sweden", "Iceland", "Liechtenstein", "Norway", "Switzerland", "United Kingdom")

EUR_dat <- EUR_dat[EUR_dat[,"location"]%in%EEAclub,]


# remove Belgium:
EUR_dat <- EUR_dat[EUR_dat[,"location"]!="Belgium",]
EUR_dat



# incorporate seroprevalence studies
zz <- getURL("https://raw.githubusercontent.com/harlanhappydog/COVID19IFR/master/serostudies_dates.csv")
phi1_all <- read.csv(text=zz)[,1:18]

# adjustement to create "aged_70_older" variable
phi1_all[,"aged_70_older"] <- 100-phi1_all[,"aged_70_younger"]
phi1_all


# Find Beta parameters to match IR CI numbers
for(k in 1:dim(phi1_all)[1]){
	IR_CI <- as.numeric(phi1_all[k,c("IR_low", "IR_high")]/100)
	ab_param <- get.beta.par(p=c(0.025,0.975),q=IR_CI, plot=FALSE, show.output=FALSE)
	phi1_all[k,"a_par"] <- ab_param[1]
	phi1_all[k,"b_par"] <- ab_param[2]}

# impute equivalent CC and T numbers (based on beta parameter values
phi1_all[,"CC"]<-round(phi1_all[,"a_par"])
phi1_all[,"T"]<-round(phi1_all[,"a_par"]+phi1_all[,"b_par"])


# what is the implied IFR:
IFR_raw <- t(apply(rbind(1: dim(phi1_all)[1]),2,function(x){
c(100*prop.test(x=phi1_all[x,"D"], n=phi1_all[x,"P"]*phi1_all[x,c("IR_high")]/100)$conf.int[1],
100*prop.test(x=phi1_all[x,"D"], n=phi1_all[x,"P"]*phi1_all[x,c("IR_low")]/100)$conf.int[2])}))

colnames(IFR_raw) <- c("IFR_low","IFR_high")

phi1_all <- cbind(phi1_all , IFR_raw)


# prepare to  merge the two datasets together:
kprime <- dim(phi1_all)[1]

EUR_dat <- cbind(k=(kprime+1):(kprime+dim(na.omit(EUR_dat))[1]), (EUR_dat))


colnames(EUR_dat)<-c(
"k", 
"Study", 
"Location",
"date",
"CC",
"T",
"total_cases_per_million",
"aged_70_older",
"hospital_beds_per_thousand", 
"P", 
"D")


# merge two datasets:
fullEUR <- rbind(phi1_all[,c("k","Study", "Location", "date","T", "CC","P","D", "aged_70_older", "hospital_beds_per_thousand")],
EUR_dat[,c("k","Study", "Location", "date","T", "CC","P","D","aged_70_older","hospital_beds_per_thousand")])
rownames(fullEUR)<-1:dim(fullEUR)[1]

# clean up :
fullEUR$CC<-round(fullEUR$CC)
fullEUR$P<-round(fullEUR$P)
fullEUR$aged_70_older <-round(fullEUR$aged_70_older)

fullEUR$iso_code <- c("CHE","DEU","LUX","HRV","CHE", as.character(fullEUR$Study[-c(1:kprime)]))

fullEUR$sero <- c(rep(1, dim(phi1_all)[1]), rep(0, dim(fullEUR)[1]-dim(phi1_all)[1]))



## add additional covariate: days_since_first_10infections

fullEUR1<-fullEUR
fullEUR1[, "days_since_first_10infections"]<-NULL
for(x in 1:dim(fullEUR1)[1]){
	
	con <- (fullEUR1$iso_code)[x]
	
caseall<-owid[owid[,"total_cases"]>=10 & 	owid[,"continent"]=="Europe",c("date","iso_code", "location", "total_cases")]

fullEUR1[fullEUR1$iso_code%in%con, "days_since_first_10infections"]<-
as.numeric(0+as.Date(fullEUR1[fullEUR1[,"iso_code"]%in%con,"date"]  ) -
 as.Date(as.character((caseall[caseall[,"iso_code"]%in%c(con),][1,"date"]))))
}


## days between first case and social distencing:

fullEUR2<-fullEUR1
d<-NULL
for(country in unique(fullEUR2$iso_code)){
days_till_lockdown1<-as.numeric(

as.Date(na.omit(owid[owid[,"iso_code"]%in%country & owid[,"stringency_index"]>=20 , c("date","location","total_cases","stringency_index")])[1,"date"]) -

as.Date(na.omit(owid[owid[,"iso_code"]%in%country & owid[,"total_cases"]>=1 , c("date","location","total_cases","stringency_index")])[1,"date"])
)
d1<-data.frame(days_till_lockdown=days_till_lockdown1, iso_code=country)
d<-rbind(d,d1)
}
d
fullEUR3<-merge(fullEUR2, d, by="iso_code", all.x=TRUE, all.y=FALSE) 

fullEUR3 <- na.omit(data.frame(fullEUR3[order(fullEUR3$k),]))


# scale all covariates:
fullEUR4 <- fullEUR3
fullEUR4[,"k"] <- 1:dim(fullEUR4)[1]

fullEUR4$scale_aged_70_older <- scale(log(fullEUR4$aged_70_older))
fullEUR4$scale_hospital_beds_per_thousand <- scale(log(fullEUR4$hospital_beds_per_thousand))
fullEUR4$scale_days_since_first_10infections <- scale(log(fullEUR4$days_since_first_10infections))
fullEUR4$scale_days_till_lockdown <- scale(log(fullEUR4$days_till_lockdown+1))

fullEUR4
fullEUR4[,"study_names"]<- factor(paste(fullEUR4[,"k"],fullEUR4[,"Location"], sep="- "), levels=c(paste(fullEUR4[,"k"],fullEUR4[,"Location"], sep="- ")))

fullEUR4[,"study_names"]
## to print:
fullEUR <- na.omit(fullEUR4)
fullEUR <- droplevels (fullEUR)



toprint <-fullEUR[,c("Location", "date", "T", "CC", "P", "D", "aged_70_older", "hospital_beds_per_thousand" , "days_since_first_10infections", "days_till_lockdown")]

print(xtable(toprint, digits=c(rep(0,7), rep(2,(dim(toprint)[2]-7)), 0) ),include.rownames=F)



#################################################
#################################################
### WRITE OUT MODELS IN STAN
#################################################
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
	
icloglogtheta ~ uniform(0, 1);
icloglogbeta ~ uniform(0, 1);

sd_sig ~ normal(0,1);
sd_tau ~ normal(0,0.1);

for (k in 1:kprime){	
	CC[k] ~ binomial(T[k], IR[k]);	}
if(kprime<K){
for (k in (kprime+1):K){	
	CC[k] ~ binomial(T[k], 1- pow((1-IR[k]), phi[k]) );	}}

	D ~ binomial(P, IR .* IFR);		
}

"


stan_mod_cov <- "


data {
int<lower=0> kprime;
int<lower=1> K;
int<lower=0> D[K];
int<lower=0> CC[K];
int<lower=0> T[K];
int<lower=0> P[K];
real<lower=0> lambda;
vector[K] scale_aged_70_older;
vector[K] scale_hospital_beds_per_thousand;
vector[K] scale_days_since_first_10infections;
vector[K] scale_days_till_lockdown;
}

parameters {
real<lower=0> gamma;
vector<lower=1, upper=1+gamma>[K] phi;
vector[K] cloglogIFR;
vector[K] cloglogIR;
real theta;
real theta1;
real theta2;
real beta;
real beta1;
real beta2;
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
theta1 ~ normal(0, 10);
theta2 ~ normal(0, 10);
beta1 ~ normal(0, 10);
beta2 ~ normal(0, 10);

gamma ~ exponential(lambda);
phi ~  uniform(1, (1+gamma));	
cloglogIFR ~ normal(theta + theta1* scale_aged_70_older + theta2* scale_hospital_beds_per_thousand, sd_tau);
cloglogIR ~ normal(beta + beta1*scale_days_since_first_10infections + beta2*scale_days_till_lockdown, sd_sig);
	
icloglogtheta ~ uniform(0, 1);
icloglogbeta ~ uniform(0, 1);

sd_sig ~ normal(0,1);
sd_tau ~ normal(0,0.1);

for (k in 1:kprime){	
	CC[k] ~ binomial(T[k], IR[k]);	}
if(kprime<K){
for (k in (kprime+1):K){	
	CC[k] ~ binomial(T[k], 1- pow((1-IR[k]), phi[k]) );	}}
	
	D ~ binomial(P, IR .* IFR);	
	
}

"

#################################################
#################################################
### FIT MODEL WITH ONLY STERO-STUDIES
#################################################
#################################################	

OB_sero <- fullEUR[fullEUR$sero==1,]

datalist_sero <- list(kprime=sum(OB_sero$sero), K=dim(OB_sero)[1],  D=round(OB_sero$D), P=round(OB_sero$P), T=round(OB_sero$T), CC=round(OB_sero$CC), lambda=0.5)


sero_simple <- stan(model_code = stan_mod_simple,
                 iter = 5000, 
                 warmup = 500,
                 data = datalist_sero,
                 chains = 4,
                 thin = 100,
                 control = list(max_treedepth = 15, adapt_delta = 0.99)
                 )


## diagnostics
rstan::traceplot(sero_simple, par=c("icloglogtheta","icloglogbeta"))


## results

100*summary(sero_simple)$summary["icloglogtheta", c("50%", "2.5%", "97.5%")]
100*summary(sero_simple)$summary["icloglogbeta", c("50%", "2.5%", "97.5%")]


## equal-tailed credible intervals
print(summary(sero_simple)$summary[c("theta", "beta", "sd_tau", "sd_sig"), c("50%", "2.5%", "97.5%")],3)



## HPD intervals:
summary_md_ci <- function(xx){
c(md=summary(sero_simple)$summary[xx,"50%"], lower=HDIofMCMC(unlist(As.mcmc.list(sero_simple, par=xx)))[1], higher=HDIofMCMC(unlist(As.mcmc.list(sero_simple, par=xx)))[2])}


data.frame(param=c("icloglogtheta","icloglogbeta","sd_tau","sd_sig"),t(apply(cbind(c("icloglogtheta","icloglogbeta","sd_tau","sd_sig")),1,function(x) summary_md_ci(x))))


##########################################
# IFR variables
IFRvars <- cbind(apply(cbind(1:dim(OB_sero)[1]), 1, function(ii) paste("IFR[",ii,"]", sep="")))

study_names <- fullEUR[ fullEUR[,"sero"]==1,"study_names"]

meta_IFR <- data.frame(study_names, 100*t(apply(IFRvars, 1, function(xx){summary_md_ci(xx)})))

meta_IFR_plus <- rbind(meta_IFR, data.frame(study_names="Model estimate", rbind(100*summary_md_ci("icloglogtheta"))))

colnames(meta_IFR_plus) <- c("Study","IFR","lower","upper")
meta_IFR_plus<-droplevels(meta_IFR_plus)

a <- ggplot(meta_IFR_plus, aes(x= Study, y=IFR,ymax=upper,ymin=lower,size=2))

b <- a+geom_pointrange(cex=c(rep(0.45, datalist_sero$K), 1), col=c(rep("tomato", datalist_sero$kprime), rep("darkolivegreen4", datalist_sero$K-datalist_sero$kprime), "black"))

#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b+coord_flip() + scale_x_discrete(limits = rev(levels(meta_IFR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()


IFR_raw <- phi1_all[, c("IFR_low",  "IFR_high")]

IFRplot_sero <- e + geom_pointrange(data=cbind(meta_IFR_plus , rbind(IFR_raw,c(NA,NA))), aes(x=Study,ymax=IFR_high,ymin=IFR_low), lwd=.82, pch="", cex=0.3, col="lightgrey", position = position_nudge(x = -0.15))+ scale_y_continuous(breaks=seq(0,2,0.25))

##########################################


##########################################
# IR variables
IRvars <- cbind(apply(cbind(1:dim(OB_sero)[1]), 1, function(ii) paste("IR[",ii,"]", sep="")))

study_names <- fullEUR[ fullEUR[,"sero"]==1,"study_names"]

meta_IR <- data.frame(study_names, 100*t(apply(IRvars, 1, function(xx){summary_md_ci(xx)})))

meta_IR_plus <- rbind(meta_IR, data.frame(study_names="Model estimate", rbind(100*summary_md_ci("icloglogbeta"))))

colnames(meta_IR_plus) <- c("Study","IR","lower","upper")
meta_IR_plus <-droplevels(meta_IR_plus)
a <- ggplot(meta_IR_plus, aes(x= Study, y=IR,ymax=upper,ymin=lower,size=2))

b <- a+geom_pointrange(cex=c(rep(0.45, datalist_sero$K), 1), col=c(rep("tomato", datalist_sero$kprime), rep("darkolivegreen4", datalist_sero$K-datalist_sero$kprime), "black"))

#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b+coord_flip() + scale_x_discrete(limits = rev(levels(meta_IR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()


IR_raw <- phi1_all[, c("IR_low",  "IR_high")]

IRplot_sero <- e + geom_pointrange(data=cbind(meta_IR_plus , rbind(IR_raw,c(NA,NA))), aes(x=Study,ymax=IR_high,ymin=IR_low), lwd=.82, pch="", cex=0.3, col="lightgrey", position = position_nudge(x = -0.15))

##########################################
IFRplot_sero_noyaxis<- IFRplot_sero + theme(axis.title.y = element_blank(),axis.text.y = element_blank())

require(gridExtra)
grid.arrange(IRplot_sero, IFRplot_sero_noyaxis, ncol=2)







#################################################
#################################################
### FIT MODEL WITH FULL DATA
#################################################
#################################################	

OB <- (fullEUR)


datalist <- list(kprime=sum(OB$sero), K=dim(OB)[1],  D=round(OB$D), P=round(OB$P), T=round(OB$T), CC=round(OB$CC), lambda=0.5, scale_aged_70_older =c(OB$scale_aged_70_older), scale_hospital_beds_per_thousand = c(OB$scale_hospital_beds_per_thousand), scale_days_since_first_10infections =c(OB$scale_days_since_first_10infections), scale_days_till_lockdown =c(OB$scale_days_till_lockdown))



MAAD <- stan(model_code = stan_mod_cov,
                 iter = 1000, 
                 warmup = 500,
                 data = datalist,
                 chains = 4,
                 thin = 10,
                 control = list(max_treedepth = 15, adapt_delta = 0.99)
                 )
                 
                 

## diagnostics
rstan::traceplot(MAAD, par=c("icloglogtheta","icloglogbeta"))
rstan::traceplot(MAAD, par=c("theta1","theta2"))
rstan::traceplot(MAAD, par=c("beta1","beta2"))
## results

100*summary(MAAD)$summary["icloglogtheta", c("50%", "2.5%", "97.5%")]
100*summary(MAAD)$summary["icloglogbeta", c("50%", "2.5%", "97.5%")]


## equal-tailed credible intervals
print(summary(MAAD)$summary[c("theta", "theta1", "theta2", "beta", "beta1",  "beta2", "sd_tau", "sd_sig"), c("50%", "2.5%", "97.5%")],3)

main_params <- c("icloglogtheta","icloglogbeta", "theta", "theta1", "theta2", "beta", "beta1",  "beta2", "sd_tau", "sd_sig")

## HPD intervals:
summary_md_ci_MAAD <- function(xx){
c(md=summary(MAAD)$summary[xx,"50%"], lower=HDIofMCMC(unlist(As.mcmc.list(MAAD, par=xx)))[1], higher=HDIofMCMC(unlist(As.mcmc.list(MAAD, par=xx)))[2])}


data.frame(param=c(main_params),t(apply(cbind(c(main_params)),1,function(x) summary_md_ci_MAAD(x))))


##########################################
# IFR variables
IFRvars <- cbind(apply(cbind(1:dim(OB)[1]), 1, function(ii) paste("IFR[",ii,"]", sep="")))

study_names <- fullEUR[ ,"study_names"]

meta_IFR <- data.frame(study_names, 100*t(apply(IFRvars, 1, function(xx){summary_md_ci_MAAD(xx)})))


meta_IFR_plus <- rbind(meta_IFR, data.frame(study_names="Model estimate", rbind(100*summary_md_ci_MAAD("icloglogtheta"))))

colnames(meta_IFR_plus) <- c("Study", "IFR","lower","upper")

a <- ggplot(meta_IFR_plus, aes(x= Study, y=IFR,ymax=upper,ymin=lower,size=2))

b <- a+geom_pointrange(cex=c(rep(0.45, datalist$K), 1), col=c(rep("tomato", datalist$kprime), rep("darkolivegreen4", datalist$K-datalist$kprime), "black"))

#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b+coord_flip() + scale_x_discrete(limits = rev(levels(meta_IFR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()

IFRplot<-e
IFRplot
##########################################


##########################################
# IR variables
IRvars <- cbind(apply(cbind(1:dim(OB)[1]), 1, function(ii) paste("IR[",ii,"]", sep="")))

study_names <- fullEUR[,"study_names"]

meta_IR <- data.frame(study_names, 100*t(apply(IRvars, 1, function(xx){summary_md_ci_MAAD(xx)})))

meta_IR_plus <- rbind(meta_IR, data.frame(study_names="Model estimate", rbind(100*summary_md_ci_MAAD("icloglogbeta"))))

colnames(meta_IR_plus) <- c("Study","IR","lower","upper")
meta_IR_plus <-droplevels(meta_IR_plus)


a <- ggplot(meta_IR_plus, aes(x= Study, y=IR,ymax=upper,ymin=lower,size=2))

b <- a+geom_pointrange(cex=c(rep(0.45, datalist$K), 1), col=c(rep("tomato", datalist$kprime), rep("darkolivegreen4", datalist$K-datalist$kprime), "black"))

#this flips the co-ordinates so your x axis becomes your y and vice versa
cc <- b+coord_flip() + scale_x_discrete(limits = rev(levels(meta_IR_plus$Study)))

#this puts in a dotted line at the point of group difference
d <- cc
e <- d + theme_bw()

IRplot <- e 
IRplot


##########################################
IFRplot_noyaxis<- IFRplot + theme(axis.title.y = element_blank(),axis.text.y = element_blank())

require(gridExtra)
grid.arrange(IRplot, IFRplot_noyaxis, ncol=2)


diff(100*HDIofMCMC(unlist(As.mcmc.list(sero_simple, par="icloglogtheta"))))
diff(100*HDIofMCMC(unlist(As.mcmc.list(MAAD, par="icloglogtheta"))))
                 