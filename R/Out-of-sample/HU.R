##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##
##  Fitting and forecasting the HU model using the 
##  demography package
##  
##  Notes: 1) here for out-of-sample
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
cou <- "CHE"      ## "SWE "or "CHE"
sex <- "Females"  ## "Females" or "Males"
name <- paste(cou,"_HMD.Rdata",sep="")
load(name)

## select ages and years
years <- 1950:2016
ages <- 0:110
m <- length(ages)
n <- length(years)
FittingData <- extract.years(FullData, years=years)

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
## select out-of-sample scenario
BT <- 10
years.fit <- 1950:(years[n]-BT)
n.fit <- length(years.fit)
years.fore <- (years.fit[n.fit]+1):2016
n.fore <- length(years.fore)
E.fit <- E[,years%in%years.fit]
Z.fit <- Z[,years%in%years.fit]
LHAZact.fit <- log(Z.fit/E.fit)
E.fore <- E[,years%in%years.fore]
Z.fore <- Z[,years%in%years.fore]
LHAZact.fore <- log(Z.fore/E.fore)
FDM_FitData <- extract.years(FittingData,years.fit)
n.boot <- 1000      ## increase for smoother PI

##-- Functional demographic model (Hyndman & Ullah 2007) -------
smus <- smooth.demogdata(FDM_FitData)
plot(FittingData)
lines(smus)
if (sex=="Males"){
  series.hu <- "male"
}else{
  series.hu <- "female"
}

## fit HU model
FDMfit <- fdm(smus,series = series.hu)
FDMfore <- forecast(FDMfit,jumpchoice="fit",h=n.fore)
if (sex=="Males"){
  ETAfitFDM <- FDMfit$male
  ETAforeFDM <- log(FDMfore$rate$male)
}else{
  ETAfitFDM <- FDMfit$female
  ETAforeFDM <- log(FDMfore$rate$female)
}

e0obs <- g0obs <- e0fitHU <- g0fitHU <- numeric(n.fit)
for (i in 1:n.fit){
  e0obs[i] <- lifetable.mx(x = ages, mx=exp(LHAZact[,i]),sex=sex)$ex[1]
  g0obs[i] <- GINI_func(ages = ages, mx=exp(LHAZact[,i]),sex=sex)
  e0fitHU[i] <- lifetable.mx(x = ages, mx=exp(ETAfitFDM[,i]),sex=sex)$ex[1]
  g0fitHU[i] <- GINI_func(ages = ages, mx=exp(ETAfitFDM[,i]),sex=sex)
}
e0foreHU <- g0foreHU <- numeric(n.fore)
for (i in 1:n.fore){
  e0foreHU[i] <- lifetable.mx(x = ages, mx=exp(ETAforeFDM[,i]),sex=sex)$ex[1]
  g0foreHU[i] <- GINI_func(ages = ages, mx=exp(ETAforeFDM[,i]),sex=sex)
}

plot(years.fit,e0obs,ylim=range(e0obs,e0foreHU),xlim=range(years,years.fore),pch=19)
lines(years.fore,e0foreHU,col=4,lwd=2)

## Forecast: BOOTSTRAP
set.seed(2019)
E0BOOT <- G0BOOT <- matrix(NA,n.fore,n.boot)
MXBOOT <- simulate(FDMfore,bootstrap = T,nsim = n.boot)
LMXBOOT <- log(MXBOOT)
jj <- 1
for (jj in 1:n.boot){
  ## for each bootstrap, compute e0, g0
  for (i in 1:n.fore){
    E0BOOT[i,jj] <- lifetable.mx(x = ages, mx=MXBOOT[,i,jj],sex=substr(sex,1,1))$ex[1]
    G0BOOT[i,jj] <- GINI_func(ages = ages, mx=MXBOOT[,i,jj],sex=substr(sex,1,1))
  }  
  cat("bootstrap",jj,"\n")
}
rm(MXBOOT)

## save results
if (BT==10){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT10y") 
}else if (BT==20){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT20y")
}else if (BT==30){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT30y")
}
name <- paste(cou,sex,"HU.Rdata",sep="_")
save.image(name)

## SOME PLOTS ## SOME PLOTS ## SOME PLOTS 
## SOME PLOTS ## SOME PLOTS ## SOME PLOTS 

## summary measures
lev.p <- 0.8
par(mar=c(4, 4, 1.5, 1.5)+0.1,mfrow=c(1,2))
plot(years.fit,e0obs,ylim=range(e0obs,E0BOOT,finite=T),xlim=range(years,years.fore),
     pch=16,xlab="Years",ylab="e(0)",main=paste("LE - ",cou,sex))
matlines(years.fore,E0BOOT,col="grey70",lwd=0.8)
lines(years.fore,e0foreHU,col=3,lwd=2)
lines(years.fore,apply(E0BOOT,1,median,na.rm=T),col=4,lwd=2)
lines(years.fore,apply(E0BOOT,1,quantile, prob=1- (1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
lines(years.fore,apply(E0BOOT,1,quantile, prob=(1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
plot(years.fit,g0obs,ylim=range(g0obs,G0BOOT,finite=T),xlim=range(years,years.fore),
     pch=16,xlab="Years",ylab="g(0)",main=paste("Gini - ",cou,sex))
matlines(years.fore,G0BOOT,col="grey70",lwd=0.8)
lines(years.fore,apply(G0BOOT,1,median,na.rm=T),col=4,lwd=2)
lines(years.fore,apply(G0BOOT,1,quantile, prob=1- (1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
lines(years.fore,apply(G0BOOT,1,quantile, prob=(1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
lines(years.fore,g0foreHU,col=3,lwd=2)
par(mfrow=c(1,1))

## END