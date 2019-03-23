

## EVALUATE BT

## cleaning the workspace
rm(list=ls(all=TRUE))

## BT scenario
BT <- 30
if (BT==10){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT10y") 
}else if (BT==20){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT20y")
}else if (BT==30){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT30y")
}

## select Population of interest
cou <- "SWE"
sex <- "Males"

## load STAD data
name <- paste(cou,sex,"3C-STAD.Rdata",sep="_")
load(name)
LMXBOOT_STAD <- STAD$LMX_STAD_fore_boot
E0BOOT_STAD <- E0BOOT
G0BOOT_STAD <- G0BOOT
var.keep <- c("LMXBOOT_STAD","E0BOOT_STAD","G0BOOT_STAD",
              "ages","n.fore","years.fore",
              "cou","sex")
rm(list=setdiff(ls(),var.keep))

## load CODA data
name <- paste(cou,sex,"CODA.Rdata",sep="_")
load(name)
var.keep <- c("LMXBOOT_STAD","E0BOOT_STAD","G0BOOT_STAD",
              "CoDa.ex1",
              "ages","n.fore","years.fore",
              "cou","sex")
rm(list=setdiff(ls(),var.keep))

## load HU data
name <- paste(cou,sex,"HU.Rdata",sep="_")
load(name)
LMXBOOT_HU <- LMXBOOT
E0BOOT_HU <- E0BOOT
G0BOOT_HU <- G0BOOT
var.keep <- c("LMXBOOT_STAD","E0BOOT_STAD","G0BOOT_STAD",
              "CoDa.ex1",
              "LMXBOOT_HU","E0BOOT_HU","G0BOOT_HU",
              "ages","n.fore","years.fore",
              "cou","sex")
rm(list=setdiff(ls(),var.keep))

## load LC data (and BT data)
name <- paste(cou,sex,"LC.Rdata",sep="_")
load(name)
LMXBOOT_LC <- LMXBOOT
E0BOOT_LC <- E0BOOT
G0BOOT_LC <- G0BOOT
## bt
MX_BT <- exp(LHAZact.fore)
e0_BT <- e0BT
g0_BT <- g0BT
var.keep <- c("LMXBOOT_STAD","E0BOOT_STAD","G0BOOT_STAD",
              "LMXBOOT_LC","E0BOOT_LC","G0BOOT_LC",
              "LMXBOOT_HU","E0BOOT_HU","G0BOOT_HU",
              "LMX_LC_cent","e0_LC_cent","g0_LC_cent",
              "ages","n.fore","years.fore","ages.lc",
              "CoDa.ex1",
              "MX_BT","e0_BT","g0_BT",
              "cou","sex")
rm(list=setdiff(ls(),var.keep))

## load forecast accuracy function
setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Functions")
source("FOREaccuracy.R")

## level
lev <- 80
lev.p <- lev/100
## STAD
e0_STAD_med <- apply(E0BOOT_STAD, 1, median, na.rm=T)
e0_STAD_up <- apply(E0BOOT_STAD, 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
e0_STAD_low <- apply(E0BOOT_STAD, 1, quantile, prob=(1-lev.p)/2,na.rm=T)
g0_STAD_med <- apply(G0BOOT_STAD, 1, median, na.rm=T)
g0_STAD_up <- apply(G0BOOT_STAD, 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
g0_STAD_low <- apply(G0BOOT_STAD, 1, quantile, prob=(1-lev.p)/2,na.rm=T)
LMX_STAD_med <- LMX_STAD_up <-  LMX_STAD_low <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_STAD_med[,i] <- apply(LMXBOOT_STAD[,i,],1,median, na.rm=T)
  LMX_STAD_up[,i] <- apply(LMXBOOT_STAD[,i,], 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
  LMX_STAD_low[,i] <- apply(LMXBOOT_STAD[,i,], 1, quantile, prob=(1-lev.p)/2,na.rm=T)
  
}

## LC
e0_LC_med <- apply(E0BOOT_LC, 1, median, na.rm=T)
e0_LC_up <- apply(E0BOOT_LC, 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
e0_LC_low <- apply(E0BOOT_LC, 1, quantile, prob=(1-lev.p)/2,na.rm=T)
g0_LC_med <- apply(G0BOOT_LC, 1, median, na.rm=T)
g0_LC_up <- apply(G0BOOT_LC, 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
g0_LC_low <- apply(G0BOOT_LC, 1, quantile, prob=(1-lev.p)/2,na.rm=T)
LMX_LC_med <- LMX_LC_up <-  LMX_LC_low <- matrix(NA,length(ages.lc),n.fore)
for (i in 1:n.fore){
  LMX_LC_med[,i] <- apply(LMXBOOT_LC[,i,],1,median, na.rm=T)
  LMX_LC_up[,i] <- apply(LMXBOOT_LC[,i,], 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
  LMX_LC_low[,i] <- apply(LMXBOOT_LC[,i,], 1, quantile, prob=(1-lev.p)/2,na.rm=T)
  
}

## HU
e0_HU_med <- apply(E0BOOT_HU, 1, median, na.rm=T)
e0_HU_up <- apply(E0BOOT_HU, 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
e0_HU_low <- apply(E0BOOT_HU, 1, quantile, prob=(1-lev.p)/2,na.rm=T)
g0_HU_med <- apply(G0BOOT_HU, 1, median, na.rm=T)
g0_HU_up <- apply(G0BOOT_HU, 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
g0_HU_low <- apply(G0BOOT_HU, 1, quantile, prob=(1-lev.p)/2,na.rm=T)
LMX_HU_med <- LMX_HU_up <-  LMX_HU_low <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_HU_med[,i] <- apply(LMXBOOT_HU[,i,],1,median, na.rm=T)
  LMX_HU_up[,i] <- apply(LMXBOOT_HU[,i,], 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
  LMX_HU_low[,i] <- apply(LMXBOOT_HU[,i,], 1, quantile, prob=(1-lev.p)/2,na.rm=T)
}


## CODA
e0_CODA_med <- CoDa.ex1$ex.med
e0_CODA_up <- CoDa.ex1$ex.up80
e0_CODA_low <- CoDa.ex1$ex.lw80
g0_CODA_med <- CoDa.ex1$gx.med
g0_CODA_up <- CoDa.ex1$gx.up80
g0_CODA_low <- CoDa.ex1$gx.lw80
MX_CODA_med <- t(CoDa.ex1$mx.med)
MX_CODA_up <-  t(CoDa.ex1$mx.up80)
MX_CODA_low <- t(CoDa.ex1$mx.lw80)


## evaluate accuracy
STAD_ACC <- ForecastAccuracyIndic(ages=ages,MXact = MX_BT,MXforeMEAN = exp(LMX_STAD_med),MXforeLOW = exp(LMX_STAD_low),MXforeUP = exp(LMX_STAD_up),
                                  LifeExpAct = e0_BT,LifeExpForeMEAN = e0_STAD_med,LifeExpForeLOW = e0_STAD_low,LifeExpForeUP = e0_STAD_up,
                                  GiniAct = 100*g0_BT,GiniForeMEAN = 100*g0_STAD_med, GiniForeLOW = 100*g0_STAD_low, GiniForeUP = 100*g0_STAD_up,
                                  alpha = lev.p,ndigits = 2)

LC_ACC <- ForecastAccuracyIndic(ages=ages.lc,MXact = MX_BT[ages%in%ages.lc,],MXforeMEAN = exp(LMX_LC_med),MXforeLOW = exp(LMX_LC_low),MXforeUP = exp(LMX_LC_up),
                                  LifeExpAct = e0_BT,LifeExpForeMEAN = e0_LC_med,LifeExpForeLOW = e0_LC_low,LifeExpForeUP = e0_LC_up,
                                  GiniAct = 100*g0_BT,GiniForeMEAN = 100*g0_LC_med,GiniForeLOW = 100*g0_LC_low,GiniForeUP = 100*g0_LC_up,
                                  alpha = lev.p,ndigits = 2)

CODA_ACC <- ForecastAccuracyIndic(ages=ages,MXact = MX_BT,MXforeMEAN = MX_CODA_med[1:length(ages),],MXforeLOW = MX_CODA_low[1:length(ages),],MXforeUP = MX_CODA_up[1:length(ages),],
                                  LifeExpAct = e0_BT,LifeExpForeMEAN = e0_CODA_med,LifeExpForeLOW = e0_CODA_low,LifeExpForeUP = e0_CODA_up,
                                  GiniAct = 100*g0_BT,GiniForeMEAN = 100*g0_CODA_med,GiniForeLOW = 100*g0_CODA_low,GiniForeUP = 100*g0_CODA_up,
                                  alpha = lev.p,ndigits = 2)

HU_ACC <- ForecastAccuracyIndic(ages=ages,MXact = MX_BT,MXforeMEAN = exp(LMX_HU_med),MXforeLOW = exp(LMX_HU_low),MXforeUP = exp(LMX_HU_up),
                                  LifeExpAct = e0_BT,LifeExpForeMEAN = e0_HU_med,LifeExpForeLOW = e0_HU_low,LifeExpForeUP = e0_HU_up,
                                  GiniAct = 100*g0_BT,GiniForeMEAN = 100*g0_HU_med,GiniForeLOW = 100*g0_HU_low,GiniForeUP = 100*g0_HU_up,
                                  alpha = lev.p,ndigits = 2)

df <- cbind(STAD_ACC$accuracy,LC_ACC$accuracy,CODA_ACC$accuracy,HU_ACC$accuracy)
colnames(df) <- c("STAD","LC","CODA","HU")
df
