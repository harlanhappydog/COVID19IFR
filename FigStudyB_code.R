############################################################
#### Bayesian adjustment for preferential testing ####
#### in estimating the COVID-19 infection fatality rate ####
#### code written by Harlan Campbell ####
#### contact harlan.campbell@stat.ubc.ca ####

#### Code to reproduce simulation study figure (Study B):

library(ggplot2)
library(gridExtra)

sim_results<-readRDS(url("https://github.com/harlanhappydog/COVID19IFR/blob/master/STUDY_B_results.rds?raw=true"))

sim_results


RR <- rbind(cbind(M=2,sim_results[,1:4]),cbind(M=3,sim_results[,6:9]),cbind(M=1,sim_results[,10:13]))
library(ggplot2, gridExtra)
RR$lambda<-as.factor(RR$lambda)
RR$gamma1<-RR$gamma+1

RR2 <-RR

g1 <- ggplot(mapping=aes(x = gamma1, y = coverage, col = lambda), data = RR[RR$M==2,]) + 
  labs(col = expression(lambda)) +  geom_line(lwd=1.1)  + geom_point(cex=2) + 
  scale_y_continuous(breaks=c(0,0.5,0.8,0.9,1), limits=c(0,1), name = expression(paste("coverage of 90% CI for ", theta))) + scale_x_continuous(trans='log10', breaks = RR$gamma1, labels = RR$gamma, name = expression(gamma))+ guides(color = FALSE)+geom_hline(yintercept = unique(RR[RR$M==3,"coverage"]), linetype = "dotted", lwd=1.1) 
g2 <- ggplot(mapping=aes(x=gamma1, y=width, col= lambda), data=RR[RR$M==2,]) + labs(col = expression(lambda)) +
  geom_line(lwd=1.1)  + geom_point(cex=2) + 
  scale_x_continuous(trans='log10', breaks=RR$gamma1, labels = RR$gamma, name=expression(gamma)) + 
  geom_abline(aes(slope = 0, intercept = unique(RR[RR$M==3,"width"]), linetype = "M3"), colour = "black", lwd = 1.1) +
  geom_abline(aes(slope = 0, intercept = 0, linetype = "M2"), colour = "black", lwd = 1.1)  + 
  geom_abline(aes(slope = 0, intercept = 0, linetype = "M1"), colour = "black", lwd = 1.1) +
  guides(colour   = FALSE, linetype = guide_legend(override.aes = list(linetype = c("dashed", "solid", "dotted")))) + 
  scale_linetype_manual(name = "Model", values = c("M1" = "dashed", "M3" = "dotted", "M2" = "solid"))

g1a <- ggplot(mapping=aes(x = gamma1, y = coverage, col = lambda), data = RR[RR$M==1,]) + 
  labs(col = expression(lambda)) +  geom_line(linetype = "dashed", lwd=1.1)  + 
  geom_point(cex=2) + scale_y_continuous(breaks=c(0,0.5,0.8,0.9,1), limits=c(0,1), name = expression(paste("coverage of 90% CI for ", theta))) + 
  scale_x_continuous(trans='log10', breaks = RR$gamma1, labels = RR$gamma, name = expression(gamma)) + 
  guides(color = FALSE)

g2a <- ggplot(mapping=aes(x=gamma1, y=width, col= lambda), data=RR[RR$M==1,]) + 
  labs(col = expression(lambda)) + geom_line(linetype = "dashed", lwd=1.1)  + 
  geom_point(cex=2) + 
  scale_x_continuous(trans='log10', breaks=RR$gamma1, labels = RR$gamma, name=expression(gamma))


grid.arrange(g1a, g2a, g1, g2, ncol=2)
