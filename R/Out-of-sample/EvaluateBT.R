##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##
##  R code to compare the forecast accuracy of different 
##  models (Table 2 in the chapter)
##  
##  Authors: Ugofilippo Basellini
##           Giancarlo Camarda
##  Last update: 24/03/2019
##
##---------------------------------------------------------##

## cleaning the workspace
rm(list=ls(all=TRUE))

## BT scenario
BT <- 10
if (BT==10){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT10y") 
}else if (BT==20){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT20y")
}else if (BT==30){
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/BT30y")
}

## select Population of interest
cou <- "CHE"      ## "CHE" or "SWE"
sex <- "Females"    ## "Females" or "Males"

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
LMXBOOT_CODA <- log(CoDa.ex1$MXBOOT)
LMXBOOT_CODA <- aperm(LMXBOOT_CODA,c(2,1,3))
LMXBOOT_CODA <- LMXBOOT_CODA[1:length(ages),,]
E0BOOT_CODA <- CoDa.ex1$E0BOOT
G0BOOT_CODA <- CoDa.ex1$G0BOOT
var.keep <- c("LMXBOOT_STAD","E0BOOT_STAD","G0BOOT_STAD",
              "LMXBOOT_CODA","E0BOOT_CODA","G0BOOT_CODA",
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
              "LMXBOOT_CODA","E0BOOT_CODA","G0BOOT_CODA",
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

## Out-of-sample (backtest) data
MX_BT <- exp(LHAZact.fore)
e0_BT <- e0BT
g0_BT <- g0BT
var.keep <- c("LMXBOOT_STAD","E0BOOT_STAD","G0BOOT_STAD",
              "LMXBOOT_LC","E0BOOT_LC","G0BOOT_LC",
              "LMXBOOT_HU","E0BOOT_HU","G0BOOT_HU",
              "LMXBOOT_CODA","E0BOOT_CODA","G0BOOT_CODA",
              "LMX_LC_cent","e0_LC_cent","g0_LC_cent",
              "ages","n.fore","years.fore","ages.lc",
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
e0_CODA_med <- apply(E0BOOT_CODA, 1, median, na.rm=T)
e0_CODA_up <- apply(E0BOOT_CODA, 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
e0_CODA_low <- apply(E0BOOT_CODA, 1, quantile, prob=(1-lev.p)/2,na.rm=T)
g0_CODA_med <- apply(G0BOOT_CODA, 1, median, na.rm=T)
g0_CODA_up <- apply(G0BOOT_CODA, 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
g0_CODA_low <- apply(G0BOOT_CODA, 1, quantile, prob=(1-lev.p)/2,na.rm=T)
LMX_CODA_med <- LMX_CODA_up <-  LMX_CODA_low <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_CODA_med[,i] <- apply(LMXBOOT_CODA[,i,],1,median, na.rm=T)
  LMX_CODA_up[,i] <- apply(LMXBOOT_CODA[,i,], 1, quantile, prob=1- (1-lev.p)/2,na.rm=T)
  LMX_CODA_low[,i] <- apply(LMXBOOT_CODA[,i,], 1, quantile, prob=(1-lev.p)/2,na.rm=T)
}

## evaluate accuracy
STAD_ACC <- ForecastAccuracyIndic(ages=ages,MXact = MX_BT,MXforeMEAN = exp(LMX_STAD_med),MXforeLOW = exp(LMX_STAD_low),MXforeUP = exp(LMX_STAD_up),
                                  LifeExpAct = e0_BT,LifeExpForeMEAN = e0_STAD_med,LifeExpForeLOW = e0_STAD_low,LifeExpForeUP = e0_STAD_up,
                                  GiniAct = 100*g0_BT,GiniForeMEAN = 100*g0_STAD_med, GiniForeLOW = 100*g0_STAD_low, GiniForeUP = 100*g0_STAD_up,
                                  alpha = lev.p,ndigits = 2)

LC_ACC <- ForecastAccuracyIndic(ages=ages.lc,MXact = MX_BT[ages%in%ages.lc,],MXforeMEAN = exp(LMX_LC_med),MXforeLOW = exp(LMX_LC_low),MXforeUP = exp(LMX_LC_up),
                                  LifeExpAct = e0_BT,LifeExpForeMEAN = e0_LC_med,LifeExpForeLOW = e0_LC_low,LifeExpForeUP = e0_LC_up,
                                  GiniAct = 100*g0_BT,GiniForeMEAN = 100*g0_LC_med,GiniForeLOW = 100*g0_LC_low,GiniForeUP = 100*g0_LC_up,
                                  alpha = lev.p,ndigits = 2)

CODA_ACC <- ForecastAccuracyIndic(ages=ages,MXact = MX_BT,MXforeMEAN = exp(LMX_CODA_med),MXforeLOW = exp(LMX_CODA_low),MXforeUP = exp(LMX_CODA_up),
                                LifeExpAct = e0_BT,LifeExpForeMEAN = e0_CODA_med,LifeExpForeLOW = e0_CODA_low,LifeExpForeUP = e0_CODA_up,
                                GiniAct = 100*g0_BT,GiniForeMEAN = 100*g0_CODA_med,GiniForeLOW = 100*g0_CODA_low,GiniForeUP = 100*g0_CODA_up,
                                alpha = lev.p,ndigits = 2)


HU_ACC <- ForecastAccuracyIndic(ages=ages,MXact = MX_BT,MXforeMEAN = exp(LMX_HU_med),MXforeLOW = exp(LMX_HU_low),MXforeUP = exp(LMX_HU_up),
                                  LifeExpAct = e0_BT,LifeExpForeMEAN = e0_HU_med,LifeExpForeLOW = e0_HU_low,LifeExpForeUP = e0_HU_up,
                                  GiniAct = 100*g0_BT,GiniForeMEAN = 100*g0_HU_med,GiniForeLOW = 100*g0_HU_low,GiniForeUP = 100*g0_HU_up,
                                  alpha = lev.p,ndigits = 2)

## compute Dawid-Sebastiani score
DSS_e0_STAD <- round(DSSfun(e0_BT,E0BOOT_STAD),2)
DSS_e0_LC <- round(DSSfun(e0_BT,E0BOOT_LC),2)
DSS_e0_CODA <- round(DSSfun(e0_BT,E0BOOT_CODA),2)
DSS_e0_HU <- round(DSSfun(e0_BT,E0BOOT_HU),2)
DSS_g0_STAD <- round(DSSfun(100*g0_BT,100*G0BOOT_STAD),2)
DSS_g0_LC <- round(DSSfun(100*g0_BT,100*G0BOOT_LC),2)
DSS_g0_CODA <- round(DSSfun(100*g0_BT,100*G0BOOT_CODA),2)
DSS_g0_HU <- round(DSSfun(100*g0_BT,100*G0BOOT_HU),2)
DSS_mx_STAD <- round(DSSmatFun(MX_BT,exp(LMXBOOT_STAD)),2)
DSS_mx_LC <- round(DSSmatFun(MX_BT[1:nrow(LMXBOOT_LC),],exp(LMXBOOT_LC)),2)
DSS_mx_CODA <- round(DSSmatFun(MX_BT,exp(LMXBOOT_CODA)),2)
DSS_mx_HU <- round(DSSmatFun(MX_BT,exp(LMXBOOT_HU)),2)


## results for the chapter
df <- cbind(STAD_ACC$accuracy,LC_ACC$accuracy,CODA_ACC$accuracy,HU_ACC$accuracy)
colnames(df) <- c("STAD","LC","CODA","HU")
df

## reshape and attach dss
df1 <- cbind(df[1:3,],df[4:6,])
rownames(df1) <- c("e0","g0","ln(mx)")
df2 <- cbind(df1,
             c(DSS_e0_STAD,DSS_g0_STAD,DSS_mx_STAD),
             c(DSS_e0_LC,DSS_g0_LC,DSS_mx_LC),
             c(DSS_e0_CODA,DSS_g0_CODA,DSS_mx_CODA),
             c(DSS_e0_HU,DSS_g0_HU,DSS_mx_HU))
colnames(df2) <- rep(c("3C-STAD","LC","CODA","HU"),3)
colnames(df2[,9:12]) <- c("STAD","LC","CODA","HU")
df2
