#####################################################################
####--------------------- Start -------- -----------------------####
#####################################################################

## cleaning the workspace
rm(list=ls(all=TRUE))

## packages
library(demography)
library(compositions)

## load CODA functions
Functions <- "~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Functions"
setwd(Functions)
source("LifeTableFUN.R")
source("CoDaFUN.R")

## load starting data
Data <- "~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data"
setwd(Data)
cou <- "SWE"      ## "SWE "or "CHE"
sex <- "Females"  ## "Females" or "Males"
name <- paste(cou,"_HMD.Rdata",sep="")
load(name)

#####################################################################
####---------------------- Information ------------------------####
#####################################################################

# Reference period
YEAR1<-1950
YEAR2<-2016
year<- c(YEAR1:YEAR2)

# Horizon forecast (until 2050)
t<-34

# Reference period + year forecast 
year.for<-c(year[-length(year)],year[length(year)]:(year[length(year)]+t))
year.pi<-c((year[length(year)]+1):(year[length(year)]+t))

#Age-range in the HMD
age1<- 0:110
#Age-range used for the forecasts
age2<- 0:120
#Age-range used to estimate the Kannisto model parameters
ar.old1<- c(80:100)+1
#Age-range used to fit the Kannisto model
ar.old2<- c(80:120)+1
#Age-range not used to fit the Kannisto model
ar.young<- c(1:80)

#Life table a,n and radix parameters
a <- c(0.06, rep(0.5, 120)) 
n <- rep(1,121)
radix<-1

# Number of simulation used for the prediction intervals (PI)
# For the uncertainty in estimating the parameters
n.error <- 10
# For the uncertainty in the extrapolated values of the time index
n.kt <- 25

#Percentile used for the PI
prob<-c(0.025, 0.1, 0.5, 0.9, 0.975)

# ARIMA order to extrapolate the time-index
# for the CoDa, with drift
order.coda<-c(0,1,1)   


#####################################################################
####-------------------------- Load data ------------------------####
#####################################################################

#######------select data for the reference period-----#########

years <- YEAR1:YEAR2
ages <- 0:110
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

Mxdta <- t(exp(LHAZact))
Exdta <- t(E)
Dxdta <- t(Z)

#####################################################################
####-------------------- Smooth and Remove 0 --------------------####
#####################################################################

#Apply the smoothing procedure
MxSmt <-fit120(Dx=Dxdta, Ex=Exdta, mx=Mxdta)

#####################################################################
####------------------- Life table calculation  -----------------####
#####################################################################

#Country-specific life tables
LTcountry <- LifeT_mx(MxSmt, radix, a ,n)

#####################################################################
####-------------------- Results --------------------------------####
#####################################################################

CoDa.ex1<-CoDa(dx=LTcountry$dx, t, sex)

setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/Fore2050")
name <- paste(cou,sex,"CODA.Rdata",sep="_")
save.image(name)

## SOME PLOTS ## SOME PLOTS ## SOME PLOTS 
## SOME PLOTS ## SOME PLOTS ## SOME PLOTS 

plot(years,LTcountry$ex[,1],xlim=range(year.for),ylim=range(CoDa.ex1$ex[,1]))
lines(year.pi,CoDa.ex1$ex.med)
lines(year.pi,CoDa.ex1$ex.lw80,col=2,lty=2)
lines(year.pi,CoDa.ex1$ex.up80,col=2,lty=2)

## END