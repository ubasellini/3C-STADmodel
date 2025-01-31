##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##
##  This R file contains the second step of the 3C-STAD model,
##  namely fitting anf forecasting component-specific 
##  age-at-death distributions with ad hoc modifications of 
##  the STAD model proposed by Basellini and Camarda (2019).
##  
##  Notes: 1) need to run the 01_FitSSE.R code first 
##         2) here for out-of-sampe
##  
##  Authors: Ugofilippo Basellini
##           Giancarlo Camarda
##  Last update: 22/03/2019
##
##---------------------------------------------------------##

## cleaning the workspace
rm(list=ls(all=TRUE))

## loading packages
library(demography)

## load functions to fit and forecast the 3C-STAD model
Functions <- "~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Functions"
setwd(Functions)
source("STAD3Cfunc.R")
source("FitSTAD3C.R")

## load starting data: observed and SSE fit 
Data <- "~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data"
setwd(Data)
cou <- "CHE"        ## "SWE "or "CHE"
sex <- "Females"    ## "Females" or "Males"
name <- paste(cou,sex,"SSEfit.Rdata",sep="_")
load(name)

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
stad.boot <- 1000      ## increase for smoother PI

## FIT AND FORE STAD
set.seed(2019)
STAD <- FitFore_STAD3C(ages=ages,years=years.fit,E=E.fit,Z=Z.fit,
                       LMX_SSE=SSEfit$LHAZhat[,years%in%years.fit],
                       xCHI=SSEfit$xCHI,LHAZ1=SSEfit$LHAZ1hat[,years%in%years.fit],
                       xSEN=SSEfit$xSEN,LHAZ2=SSEfit$LHAZ2hat[,years%in%years.fit],
                       xADU=SSEfit$xADU,LHAZ3=SSEfit$LHAZ3hat[,years%in%years.fit],
                       fore.year=years.fore[n.fore],
                       stad.boot=stad.boot,
                       PLOT=F,BOOT=T,printFIT=T)

E0BOOT <- STAD$e0stad_fore_boot
G0BOOT <- STAD$g0stad_fore_boot

if (BT==10){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT10y") 
}else if (BT==20){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT20y")
}else if (BT==30){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT30y")
}
name <- paste(cou,sex,"3C-STAD.Rdata",sep="_")
save.image(name)

## SOME PLOTS ## SOME PLOTS ## SOME PLOTS 
## SOME PLOTS ## SOME PLOTS ## SOME PLOTS 

# load(name)
## summary measures
lev.p <- 0.8
par(mar=c(4, 4, 1.5, 1.5)+0.1,mfrow=c(1,2))
plot(years.fit,SSEfit$e0act[years%in%years.fit],ylim=range(SSEfit$e0act,E0BOOT,finite=T),xlim=range(years,years.fore),
     pch=16,xlab="Years",ylab="e(0)",main=paste("LE - ",cou,sex))
matlines(years.fore,E0BOOT,col="grey70",lwd=0.8)
lines(years.fore,apply(E0BOOT,1,median,na.rm=T),col=4,lwd=2)
lines(years.fore,apply(E0BOOT,1,quantile, prob=1- (1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
lines(years.fore,apply(E0BOOT,1,quantile, prob=(1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
points(years.fore,SSEfit$e0act[years%in%years.fore])
legend("topleft",c("observed","STAD"),pch=c(16,4),col=c(1,4),
       bty="n",lwd=2,cex=1.3,lty=NA)
plot(years.fit,SSEfit$g0act[years%in%years.fit],ylim=range(SSEfit$g0act,G0BOOT,finite=T),xlim=range(years,years.fore),
     pch=16,xlab="Years",ylab="g(0)",main=paste("Gini - ",cou,sex))
matlines(years.fore,G0BOOT,col="grey70",lwd=0.8)
lines(years.fore,apply(G0BOOT,1,median,na.rm=T),col=4,lwd=2)
lines(years.fore,apply(G0BOOT,1,quantile, prob=1- (1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
lines(years.fore,apply(G0BOOT,1,quantile, prob=(1-lev.p)/2,na.rm=T),col=4,lwd=2,lty=2)
lines(years.fore,STAD$g0stad_fore,col=3,lwd=2)
points(years.fore,SSEfit$g0act[years%in%years.fore])
par(mfrow=c(1,1))

LMX_STAD_med <- LMX_STAD_up <-  LMX_STAD_low <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_STAD_med[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,median, na.rm=T)
  LMX_STAD_up[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,], 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
  LMX_STAD_low[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,], 1, quantile, prob=(1-lev.p)/2,na.rm=T)
}

## PARS ## PARS ## 
par(mfrow=c(1,3))
plot(years.fit,STAD$SENmass,ylim=range(1,STAD$SENmass),xlim=range(years,years.fore))
abline(h=1)
matlines(years.fore,STAD$SENmassBOOT,col="grey80",lty=1)
lines(years.fore,apply(STAD$SENmassBOOT,1,median,col=4))
plot(years.fit,STAD$CHImass,ylim=range(0,STAD$CHImass),xlim=range(years,years.fore))
abline(h=0)
matlines(years.fore,STAD$CHImassBOOT,col="grey80",lty=1)
lines(years.fore,apply(STAD$CHImassBOOT,1,median,col=4))
plot(years.fit,STAD$ADUmass,ylim=range(0,STAD$ADUmass,STAD$ADUmassBOOT),xlim=range(years,years.fore))
abline(h=0)
matlines(years.fore,STAD$ADUmassBOOT,col="grey80",lty=1)
lines(years.fore,apply(STAD$ADUmassBOOT,1,median,col=4))
par(mfrow=c(1,1))

## SENESCENT: BL, BU
plot(years.fit,STAD$bsLO_sen,ylim=range(STAD$bsLO_sen,STAD$BLsenBOOT),xlim=range(years,years.fore))
matlines(years.fore,STAD$BLsenBOOT,col="grey80",lty=1)
lines(years.fore,STAD$BLsen_fore)
lines(years.fore,apply(STAD$BLsenBOOT,1,median),col=4)
plot(years.fit,STAD$bsUP_sen,ylim=range(STAD$bsUP_sen,STAD$BUsenBOOT),xlim=range(years,years.fore))
matlines(years.fore,STAD$BUsenBOOT,col="grey80",lty=1)
lines(years.fore,STAD$BUsen_fore)
lines(years.fore,apply(STAD$BUsenBOOT,1,median),col=4)


## ADULTHOOD
plot(years.fit,STAD$bs_adu,ylim=range(STAD$bs_adu,STAD$BaduBOOT),xlim=range(years,years.fore))
matlines(years.fore,STAD$BaduBOOT,col="grey80",lty=1)
lines(years.fore,STAD$Badu_fore)
lines(years.fore,apply(STAD$BaduBOOT,1,median),col=4)

## END