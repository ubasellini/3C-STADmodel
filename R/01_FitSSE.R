##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##
##  This R file contains the first step of the 3C-STAD model,
##  namely decomposing the mortality pattern via the SSE
##  model proposed by Camarda et al. (2016).
##  
##  Authors: Ugofilippo Basellini
##           Giancarlo Camarda
##  Last update: 20/03/2019
##
##---------------------------------------------------------##

## cleaning the workspace
rm(list=ls(all=TRUE))

## loading packages
library(demography)
library(colorspace)

## load function to fit the SSE model
Functions <- "~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Functions"
setwd(Functions)
source("FitSSE2D.R")

## load country-specific data
setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data")
cou <- "SWE"      ## "CHE" or "SWE"
sex <- "Males"    ## "Females" or "Males"

## load country-specific data
name <- paste(cou,"_HMD.Rdata",sep="")
load(name)

## select ages, years and subset dataset
ages <- 0:110
years <- 1950:2016
m <- length(ages)
n <- length(years)
coly <- rainbow_hcl(n)
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

## SSE DECOMP
Start_Time <- Sys.time()
SSEfit <- Fit_SSE2D(ages0 = ages,years = years,Z0=Z,E0=E,sex=sex)
End_Time <- Sys.time()
(Tot_Time <- round(End_Time-Start_Time))

## save data
name <- paste(cou,sex,"SSEfit.Rdata",sep="_")
setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data")
save.image(name)

## SOME PLOTS ## SOME PLOTS ## SOME PLOTS 
## SOME PLOTS ## SOME PLOTS ## SOME PLOTS 

## summary measures
par(mfrow=c(1,2),mar=c(3,3,1,1))
plot(years,SSEfit$e0act,main="e0",pch=16)
points(years,SSEfit$e0hat,col=2,lwd=2,pch=4)
legend("topleft",c("Obs","SSE 2D"),pch=c(16,4),col=1:2,bty="n",cex=1.5,lwd=c(1,2),lty=NA)
plot(years,SSEfit$g0act,main="Gini",pch=16)
points(years,SSEfit$g0hat,col=2,lwd=2,pch=4)
par(mfrow=c(1,1))

## log-mortality rates
i <- 1
for (i in 1:n){
  par(mar=c(4, 4, 1.5, 1.5) +0.1)
  plot(ages,LHAZact[,i],t="p",col=coly[i],pch=1,main=paste(cou,sex,"-",years[i]),ylim=c(-12,0),
       xlab="Ages",ylab="Log(m(x))")
  lines(ages,SSEfit$LHAZhat[,i],col=1,lwd=2)
  lines(SSEfit$xCHI,SSEfit$LHAZ1hat[,i],col=2,lty=2,lwd=2)
  lines(SSEfit$xSEN,SSEfit$LHAZ2hat[,i],col=3,lty=2,lwd=2)
  lines(SSEfit$xADU,SSEfit$LHAZ3hat[,i],col=4,lty=2,lwd=2)
  Sys.sleep(0.2)
}

## END