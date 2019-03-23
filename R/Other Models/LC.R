##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##
##  Fitting and forecasting the LC model using the 
##  demography package
##  
##  Authors: Ugofilippo Basellini
##           Giancarlo Camarda
##  Last update: 22/03/2019
##
##---------------------------------------------------------##

## cleaning the workspace
rm(list=ls(all=TRUE))

## loading package
library(demography)

## load block diagonal function
Functions <- "~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Functions"
setwd(Functions)
source("LifeTableFUN.R")

## load starting data
Data <- "~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data"
setwd(Data)
cou <- "SWE"      ## "SWE "or "CHE"
sex <- "Females"  ## "Females" or "Males"
name <- paste(cou,"_HMD.Rdata",sep="")
load(name)

## select ages and years
years <- 1950:2016
ages <- 0:110
m <- length(ages)
n <- length(years)
FittingData <- extract.years(FullData, years=years)
LC_FitData <- FittingData
if (sex=="Males"){
  E <- FittingData$pop$male
  MX <- FittingData$rate$male
  Z <- E*MX
  LHAZact <- log(Z/E)
}else if (sex=="Females"){
  E <- FittingData$pop$female
  MX <- FittingData$rate$female
  Z <- E*MX
  LHAZact <- log(Z/E)
}

## forecast years
years.fore <- (years[n]+1):2050
n.fore <- length(years.fore)
lc.boot <- 250

##-- ORIGINAL LEE-CARTER -------
if (sex=="Males"){
  series.lc <- "male"
}else{
  series.lc <- "female"
}

iter <- 0
for (iter in 0:m){
  ages.lc <- ages[1]:ages[m-iter]
  m.lc <- length(ages.lc)
  ## tro to fit
  LCfit <- try(lca(LC_FitData, adjust="dt",years = years,
             ages = ages.lc,interpolate = T,series = series.lc,
             max.age = ages.lc[m.lc]),silent=T)
  if (class(LCfit)!="try-error") break
}

## LC pars
AlphaLC <- LCfit$ax
BetaLC <- LCfit$bx
KappaLC <- LCfit$kt

par(mfrow=c(1,3))
plot(ages.lc,AlphaLC,main = expression(alpha[x]))
plot(ages.lc,BetaLC,main = expression(beta[x]))
plot(years,KappaLC,main = expression(kappa[t]))
par(mfrow=c(1,1))

## central forecast
Kappats <- ts(c(KappaLC), start = years[1])
modKappa <- Arima(Kappats, order=c(0,1,0), include.drift=TRUE)
predK <- forecast(modKappa,h=n.fore)

## compute log-mortality LC
LCetaFUN <- function(alpha,beta,kappa){
  One <- matrix(1, nrow=length(kappa), ncol=1)
  p1 <- alpha %*% t(One)
  p2 <- beta %*% t(kappa)
  eta <- p1+p2
  return(eta)
}
ETAhatLC <- LCetaFUN(AlphaLC,BetaLC,KappaLC)
ETAforeLC <- LCetaFUN(AlphaLC,BetaLC,predK$mean)

## summary measures
e0obs <- g0obs <- e0hatLC <- g0hatLC <- numeric(n)
for (i in 1:n){
  e0obs[i] <- lifetable.mx(x = ages, mx=exp(LHAZact[,i]),sex=sex)$ex[1]
  g0obs[i] <- GINI_func(ages = ages, mx=exp(LHAZact[,i]),sex=sex)
  e0hatLC[i] <- lifetable.mx(x = ages.lc, mx=exp(ETAhatLC[,i]),sex=sex)$ex[1]
  g0hatLC[i] <- GINI_func(ages = ages.lc, mx=exp(ETAhatLC[,i]),sex=sex)
}
e0foreLC <- g0foreLC <- numeric(n.fore)
for (i in 1:n.fore){
  e0foreLC[i] <- lifetable.mx(x = ages.lc, mx=exp(ETAforeLC[,i]),sex=sex)$ex[1]
  g0foreLC[i] <- GINI_func(ages = ages.lc, mx=exp(ETAforeLC[,i]),sex=sex)
}
plot(years,e0obs,pch=19,
     xlim=range(years,years.fore),ylim=range(e0obs,e0hatLC,e0foreLC))
points(years,e0hatLC,col=4,pch=4,lwd=2)
lines(years.fore,e0foreLC,col=4,lwd=2)

## Forecast: BOOTSTRAP
set.seed(2019)
LMXBOOT <- array(NA,dim=c(m.lc,n.fore,lc.boot))
E0BOOT <- G0BOOT <- matrix(NA,n.fore,lc.boot)

## Bootstrap of kappa_t 
jj <- 1
for(jj in 1:lc.boot){
  ## generate simulation with bootsrapping
  k.sim <- simulate(modKappa, nsim=n.fore,
                    future=TRUE, bootstrap=TRUE)
  LMXBOOT[,,jj] <- LCetaFUN(alpha = AlphaLC,beta=BetaLC,
                            kappa = k.sim)
  for (i in 1:n.fore){
    E0BOOT[i,jj] <- lifetable.mx(x = ages.lc, mx=exp(LMXBOOT[,i,jj]),sex=substr(sex,1,1))$ex[1]
    G0BOOT[i,jj] <- GINI_func(ages = ages.lc, mx=exp(LMXBOOT[,i,jj]),sex=substr(sex,1,1))
  }  
  cat("bootstrap",jj,"\n")
}

setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/Fore2050")
name <- paste(cou,sex,"LC.Rdata",sep="_")
save.image(name)

## SOME PLOTS ## SOME PLOTS ## SOME PLOTS 
## SOME PLOTS ## SOME PLOTS ## SOME PLOTS 

## summary measures
lev.p <- 0.8
par(mar=c(4, 4, 1.5, 1.5)+0.1,mfrow=c(1,2))
plot(years,e0obs,ylim=range(e0obs,E0BOOT,finite=T),xlim=range(years,years.fore),
     pch=16,xlab="Years",ylab="e(0)",main=paste("LE - ",cou,sex))
points(years,e0hatLC,col=4,pch=4,lwd=2)
matlines(years.fore,E0BOOT,col="grey70",lwd=0.8)
lines(years.fore,e0foreLC,col=3,lwd=2)
lines(years.fore,apply(E0BOOT,1,median,na.rm=T),col=4,lwd=2)
lines(years.fore,apply(E0BOOT,1,quantile, prob=1- (1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
lines(years.fore,apply(E0BOOT,1,quantile, prob=(1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
plot(years,g0obs,ylim=range(g0obs,G0BOOT,finite=T),xlim=range(years,years.fore),
     pch=16,xlab="Years",ylab="g(0)",main=paste("Gini - ",cou,sex))
points(years,g0hatLC,col=4,pch=4,lwd=2)
matlines(years.fore,G0BOOT,col="grey70",lwd=0.8)
lines(years.fore,apply(G0BOOT,1,median,na.rm=T),col=4,lwd=2)
lines(years.fore,apply(G0BOOT,1,quantile, prob=1- (1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
lines(years.fore,apply(G0BOOT,1,quantile, prob=(1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
lines(years.fore,g0foreLC,col=3,lwd=2)
par(mfrow=c(1,1))

LMX_LC_med <- LMX_LC_up <-  LMX_LC_low <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_LC_med[,i] <- apply(LMXBOOT[,i,],1,median, na.rm=T)
  LMX_LC_up[,i] <- apply(LMXBOOT[,i,], 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
  LMX_LC_low[,i] <- apply(LMXBOOT[,i,], 1, quantile, prob=(1-lev.p)/2,na.rm=T)
}

matplot(ages,LHAZact,t="l",lty=1,col="grey80",ylim=c(-12,0))
lines(ages,LMX_LC_med[,n.fore],col=4,lwd=2)
lines(ages,LMX_LC_up[,n.fore],col=4,lwd=2,lty=2)
lines(ages,LMX_LC_low[,n.fore],col=4,lwd=2,lty=2)
lines(ages,ETAforeLC[,n.fore],col=2,lwd=2)

## END