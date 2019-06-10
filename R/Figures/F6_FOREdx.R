##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##
##  R code to generate Figure 6 of the chapter 
##  Forecast age-at-death distribution
##  
##  Authors: Ugofilippo Basellini
##           Giancarlo Camarda
##  Last update: 24/03/2019
##
##---------------------------------------------------------##

## cleaning the workspace
rm(list=ls(all=TRUE))

##---- SWEDISH MALES  ------
cou <- "SWE"
sex <- "Males"

## load STAD data
DataRes <- "~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Data/Fore2050"
setwd(DataRes)
name <- paste(cou,sex,"3C-STAD.Rdata",sep="_")
load(name)
lev.p <- 0.8
LMX_STAD_SWEmales <- LMX_STAD_SWEmales_up <-  LMX_STAD_SWEmales_low <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_STAD_SWEmales[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,median, na.rm=T)
  LMX_STAD_SWEmales_up[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,quantile, prob=1- (1-lev.p)/2,na.rm=T)
  LMX_STAD_SWEmales_low[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,quantile, prob=(1-lev.p)/2,na.rm=T)
}
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load CODA data
name <- paste(cou,sex,"CODA.Rdata",sep="_")
load(name)
LMX_CODA_SWEmales <- log(t(CoDa.ex1$mx.med))
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load HU data
name <- paste(cou,sex,"HU.Rdata",sep="_")
load(name)
LMX_HU_SWEmales <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_HU_SWEmales[,i] <- apply(LMXBOOT[,i,],1,median, na.rm=T)
}
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load LC data (and BT data)
name <- paste(cou,sex,"LC.Rdata",sep="_")
load(name)
ages.lc_SWEmales <- ages.lc
LMX_LC_SWEmales <- matrix(NA,length(ages.lc_SWEmales),n.fore)
for (i in 1:n.fore){
  LMX_LC_SWEmales[,i] <- apply(LMXBOOT[,i,],1,median, na.rm=T)
}
LMXobs_SWEmales <- LHAZact
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))


##---- SWEDISH FEMALES  ------
cou <- "SWE"
sex <- "Females"

## load STAD data
name <- paste(cou,sex,"3C-STAD.Rdata",sep="_")
load(name)
LMX_STAD_SWEfemales <- LMX_STAD_SWEfemales_up <-  LMX_STAD_SWEfemales_low <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_STAD_SWEfemales[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,median, na.rm=T)
  LMX_STAD_SWEfemales_up[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,quantile, prob=1- (1-lev.p)/2,na.rm=T)
  LMX_STAD_SWEfemales_low[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,quantile, prob=(1-lev.p)/2,na.rm=T)
}
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load CODA data
name <- paste(cou,sex,"CODA.Rdata",sep="_")
load(name)
LMX_CODA_SWEfemales <- log(t(CoDa.ex1$mx.med))
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low","LMX_CODA_SWEfemales",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load HU data
name <- paste(cou,sex,"HU.Rdata",sep="_")
load(name)
LMX_HU_SWEfemales <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_HU_SWEfemales[,i] <- apply(LMXBOOT[,i,],1,median, na.rm=T)
}

var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low","LMX_CODA_SWEfemales",
              "LMX_HU_SWEfemales",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load LC data (and BT data)
name <- paste(cou,sex,"LC.Rdata",sep="_")
load(name)
ages.lc_SWEfemales <- ages.lc
LMX_LC_SWEfemales <- matrix(NA,length(ages.lc_SWEfemales),n.fore)
for (i in 1:n.fore){
  LMX_LC_SWEfemales[,i] <- apply(LMXBOOT[,i,],1,median, na.rm=T)
}
LMXobs_SWEfemales <- LHAZact
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low","LMX_CODA_SWEfemales",
              "LMX_HU_SWEfemales","LMXobs_SWEfemales",
              "LMX_LC_SWEfemales","ages.lc_SWEfemales",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))


##---- CHE MALES  ------
cou <- "CHE"
sex <- "Males"

## load STAD data
name <- paste(cou,sex,"3C-STAD.Rdata",sep="_")
load(name)
LMX_STAD_CHEmales <- LMX_STAD_CHEmales_up <- LMX_STAD_CHEmales_low <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_STAD_CHEmales[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,median, na.rm=T)
  LMX_STAD_CHEmales_up[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,quantile, prob=1- (1-lev.p)/2,na.rm=T)
  LMX_STAD_CHEmales_low[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,quantile, prob=(1-lev.p)/2,na.rm=T)
}
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low","LMX_CODA_SWEfemales",
              "LMX_HU_SWEfemales","LMXobs_SWEfemales",
              "LMX_LC_SWEfemales","ages.lc_SWEfemales",
              "LMX_STAD_CHEmales","LMX_STAD_CHEmales_up","LMX_STAD_CHEmales_low",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load CODA data
name <- paste(cou,sex,"CODA.Rdata",sep="_")
load(name)
LMX_CODA_CHEmales <- log(t(CoDa.ex1$mx.med))
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low","LMX_CODA_SWEfemales",
              "LMX_HU_SWEfemales","LMXobs_SWEfemales",
              "LMX_LC_SWEfemales","ages.lc_SWEfemales",
              "LMX_STAD_CHEmales","LMX_STAD_CHEmales_up","LMX_STAD_CHEmales_low",
              "LMX_CODA_CHEmales",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load HU data
name <- paste(cou,sex,"HU.Rdata",sep="_")
load(name)
LMX_HU_CHEmales <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_HU_CHEmales[,i] <- apply(LMXBOOT[,i,],1,median, na.rm=T)
}
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low","LMX_CODA_SWEfemales",
              "LMX_HU_SWEfemales","LMXobs_SWEfemales",
              "LMX_LC_SWEfemales","ages.lc_SWEfemales",
              "LMX_STAD_CHEmales","LMX_STAD_CHEmales_up","LMX_STAD_CHEmales_low",
              "LMX_CODA_CHEmales","LMX_HU_CHEmales",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load LC data (and BT data)
name <- paste(cou,sex,"LC.Rdata",sep="_")
load(name)
ages.lc_CHEmales <- ages.lc
LMX_LC_CHEmales <- matrix(NA,length(ages.lc_CHEmales),n.fore)
for (i in 1:n.fore){
  LMX_LC_CHEmales[,i] <- apply(LMXBOOT[,i,],1,median, na.rm=T)
}
LMXobs_CHEmales <- LHAZact
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low","LMX_CODA_SWEfemales",
              "LMX_HU_SWEfemales","LMXobs_SWEfemales",
              "LMX_LC_SWEfemales","ages.lc_SWEfemales",
              "LMX_STAD_CHEmales","LMX_STAD_CHEmales_up","LMX_STAD_CHEmales_low",
              "LMX_CODA_CHEmales","LMX_HU_CHEmales","LMX_LC_CHEmales",
              "LMXobs_CHEmales","ages.lc_CHEmales",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

##---- CHE FEMALES  ------
cou <- "CHE"
sex <- "Females"

## load STAD data
name <- paste(cou,sex,"3C-STAD.Rdata",sep="_")
load(name)
LMX_STAD_CHEfemales <- LMX_STAD_CHEfemales_up <- LMX_STAD_CHEfemales_low <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_STAD_CHEfemales[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,median, na.rm=T)
  LMX_STAD_CHEfemales_up[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,quantile, prob=1- (1-lev.p)/2,na.rm=T)
  LMX_STAD_CHEfemales_low[,i] <- apply(STAD$LMX_STAD_fore_boot[,i,],1,quantile, prob=(1-lev.p)/2,na.rm=T)
  
}
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low",
              "LMX_CODA_SWEfemales",
              "LMX_HU_SWEfemales","LMXobs_SWEfemales",
              "LMX_LC_SWEfemales","ages.lc_SWEfemales",
              "LMX_STAD_CHEmales","LMX_STAD_CHEmales_up","LMX_STAD_CHEmales_low",
              "LMX_CODA_CHEmales","LMX_HU_CHEmales","LMX_LC_CHEmales",
              "LMXobs_CHEmales","ages.lc_CHEmales",
              "LMX_STAD_CHEfemales","LMX_STAD_CHEfemales_up","LMX_STAD_CHEfemales_low",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load CODA data
name <- paste(cou,sex,"CODA.Rdata",sep="_")
load(name)
LMX_CODA_CHEfemales <- log(t(CoDa.ex1$mx.med))
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low","LMX_CODA_SWEfemales",
              "LMX_HU_SWEfemales","LMXobs_SWEfemales",
              "LMX_LC_SWEfemales","ages.lc_SWEfemales",
              "LMX_STAD_CHEmales","LMX_STAD_CHEmales_up","LMX_STAD_CHEmales_low",
              "LMX_CODA_CHEmales",
              "LMX_HU_CHEmales","LMX_LC_CHEmales",
              "LMXobs_CHEmales","ages.lc_CHEmales",
              "LMX_STAD_CHEfemales","LMX_STAD_CHEfemales_up","LMX_STAD_CHEfemales_low",
              "LMX_CODA_CHEfemales",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load HU data
name <- paste(cou,sex,"HU.Rdata",sep="_")
load(name)
LMX_HU_CHEfemales <- matrix(NA,length(ages),n.fore)
for (i in 1:n.fore){
  LMX_HU_CHEfemales[,i] <- apply(LMXBOOT[,i,],1,median, na.rm=T)
}
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low","LMX_CODA_SWEfemales",
              "LMX_HU_SWEfemales","LMXobs_SWEfemales",
              "LMX_LC_SWEfemales","ages.lc_SWEfemales",
              "LMX_STAD_CHEmales","LMX_STAD_CHEmales_up","LMX_STAD_CHEmales_low",
              "LMX_CODA_CHEmales","LMX_HU_CHEmales","LMX_LC_CHEmales",
              "LMXobs_CHEmales","ages.lc_CHEmales",
              "LMX_STAD_CHEfemales","LMX_STAD_CHEfemales_up","LMX_STAD_CHEfemales_low",
              "LMX_CODA_CHEfemales",
              "LMX_HU_CHEfemales",
              "cou","sex","lev.p")
rm(list=setdiff(ls(),var.keep))

## load LC data (and BT data)
name <- paste(cou,sex,"LC.Rdata",sep="_")
load(name)
ages.lc_CHEfemales <- ages.lc
LMX_LC_CHEfemales <- matrix(NA,length(ages.lc_CHEfemales),n.fore)
for (i in 1:n.fore){
  LMX_LC_CHEfemales[,i] <- apply(LMXBOOT[,i,],1,median, na.rm=T)
}
LMXobs_CHEfemales <- LHAZact
var.keep <- c("LMX_STAD_SWEmales","LMX_STAD_SWEmales_up","LMX_STAD_SWEmales_low",
              "LMX_CODA_SWEmales",
              "LMX_HU_SWEmales","LMX_LC_SWEmales",
              "LMXobs_SWEmales","ages.lc_SWEmales",
              "LMX_STAD_SWEfemales","LMX_STAD_SWEfemales_up","LMX_STAD_SWEfemales_low","LMX_CODA_SWEfemales",
              "LMX_HU_SWEfemales","LMXobs_SWEfemales",
              "LMX_LC_SWEfemales","ages.lc_SWEfemales",
              "LMX_STAD_CHEmales","LMX_STAD_CHEmales_up","LMX_STAD_CHEmales_low",
              "LMX_CODA_CHEmales","LMX_HU_CHEmales","LMX_LC_CHEmales",
              "LMXobs_CHEmales","ages.lc_CHEmales",
              "LMX_STAD_CHEfemales","LMX_STAD_CHEfemales_up","LMX_STAD_CHEfemales_low",
              "LMX_CODA_CHEfemales",
              "LMX_HU_CHEfemales","LMX_LC_CHEfemales",
              "ages.lc_CHEfemales","LMXobs_CHEfemales","m",
              "cou","sex","lev.p","years","years.fore","n","n.fore","ages",
              "lifetable.mx")
rm(list=setdiff(ls(),var.keep))


##---- PLOTTING -----
library("RColorBrewer")
display.brewer.pal(n=8,"Blues")
col.obs <- c("grey20","grey40")
ylimdxf_SWE <- range(lifetable.mx(x=ages,mx=exp(LMXobs_SWEfemales[,n]),sex="F")$dx)
ylimLMXm_SWE <- exp(range(LMX_HU_SWEmales,LMX_LC_SWEmales,LMX_STAD_SWEmales,LMX_CODA_SWEmales,finite=T))
ylimLMXf_CHE <- exp(range(LMX_HU_CHEfemales,LMX_STAD_CHEfemales,LMX_CODA_CHEfemales,finite=T))
ylimLMXm_CHE <- exp(range(LMX_HU_CHEmales,LMX_STAD_CHEmales,LMX_CODA_CHEmales,finite=T))
ylimdx <- range(0,6500)
col.mod <- c(brewer.pal(n=8,"Blues")[7],
             brewer.pal(n=8,"Greens")[6],
             brewer.pal(n=8,"Oranges")[6],
             brewer.pal(n=8,"Purples")[8])
cex.main <- 1.3
cex.x.lab <- 1.75
cex.y.lab <- 1.5
cex.x.axis <- 1.25
cex.y.axis <- 1.1
cex.leg <- 1.3

setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/Chapter/Figures")
pdf("F6.pdf",width = 10, height = 8)
op <- par(mfrow = c(2,2),
          oma = c(2,2,0,0) + 0.1,
          mar = c(1,1.5,1,1) + 0.1)
## bottom, left, top, right
## LMX - SWE FEM
plot(ages,lifetable.mx(x=ages,mx=exp(LMXobs_SWEfemales[,n]),sex="F")$dx,t="n",main="Sweden - Females",ylim=ylimdx,
     xlab="",ylab="",axes=F,cex.main=cex.main)
axis(1,at=seq(0,100,20),labels = rep("",6))
axis(2,las=2,cex.axis=cex.y.axis,
     at=seq(0,6000,1000),labels = seq(0,6,1))
grid();box()
mtext("x 10^3", 3, cex=0.6, at=0)
mtext(text = expression(d[x]),side=2,
      line = 1.75,cex=cex.y.lab,las=1)
points(ages,lifetable.mx(x=ages,mx=exp(LMXobs_SWEfemales[,n]),sex="F")$dx,pch=1,col=col.obs[1])
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_STAD_SWEfemales[,n.fore]),sex="F")$dx,col=col.mod[1],lwd=2)
lines(ages.lc_SWEfemales,lifetable.mx(x=ages.lc_SWEfemales,mx=exp(LMX_LC_SWEfemales[,n.fore]),sex="F")$dx,col=col.mod[2],lwd=2,lty=2)
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_CODA_SWEfemales[1:m,n.fore]),sex="F")$dx,col=col.mod[3],lwd=2,lty=4)
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_HU_SWEfemales[,n.fore]),sex="F")$dx,col=col.mod[4],lwd=2,lty=5)
legend("topleft",c("Obs, 2016","3C-STAD","LC","CODA","HU"),col=c(col.obs[1],col.mod),lty=c(NA,1,2,4,5),bty="n",
       cex=cex.leg,lwd=2,pch=c(1,NA,NA,NA,NA))

## LMX - SWE MAL
plot(ages,lifetable.mx(x=ages,mx=exp(LMXobs_SWEmales[,n]),sex="F")$dx,t="n",main="Sweden - Males",ylim=ylimdx,
     xlab="",ylab="",axes=F,cex.main=cex.main)
axis(1,at=seq(0,100,20),labels = rep("",6))
axis(2,at=seq(0,6000,1000),labels = rep("",7))
grid();box()
# mtext(text = expression(log(m[x])),side=2,
#       line = 1.75,cex=cex.y.lab,las=3)
points(ages,lifetable.mx(x=ages,mx=exp(LMXobs_SWEmales[,n]),sex="F")$dx,pch=1,col=col.obs[1])
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_STAD_SWEmales[,n.fore]),sex="F")$dx,col=col.mod[1],lwd=2)
lines(ages.lc_SWEmales,lifetable.mx(x=ages.lc_SWEmales,mx=exp(LMX_LC_SWEmales[,n.fore]),sex="F")$dx,col=col.mod[2],lwd=2,lty=2)
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_CODA_SWEmales[1:m,n.fore]),sex="F")$dx,col=col.mod[3],lwd=2,lty=4)
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_HU_SWEmales[,n.fore]),sex="F")$dx,col=col.mod[4],lwd=2,lty=5)


## LMX - CHE FEM
plot(ages,lifetable.mx(x=ages,mx=exp(LMXobs_CHEfemales[,n]),sex="F")$dx,t="n",main="Switzerland - Females",ylim=ylimdx,
     xlab="",ylab="",axes=F,cex.main=cex.main)
axis(1)
axis(2,las=2,cex.axis=cex.y.axis,
     at=seq(0,6000,1000),labels = seq(0,6,1))
grid();box()
mtext("x 10^3", 3, cex=0.6, at=0)
mtext(text = expression(d[x]),side=2,
      line = 1.75,cex=cex.y.lab,las=1)
points(ages,lifetable.mx(x=ages,mx=exp(LMXobs_CHEfemales[,n]),sex="F")$dx,pch=1,col=col.obs[1])
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_STAD_CHEfemales[,n.fore]),sex="F")$dx,col=col.mod[1],lwd=2)
lines(ages.lc_CHEfemales,lifetable.mx(x=ages.lc_CHEfemales,mx=exp(LMX_LC_CHEfemales[,n.fore]),sex="F")$dx,col=col.mod[2],lwd=2,lty=2)
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_CODA_CHEfemales[1:m,n.fore]),sex="F")$dx,col=col.mod[3],lwd=2,lty=4)
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_HU_CHEfemales[,n.fore]),sex="F")$dx,col=col.mod[4],lwd=2,lty=5)

## LMX - CHE FEM
plot(ages,lifetable.mx(x=ages,mx=exp(LMXobs_CHEmales[,n]),sex="F")$dx,t="n",main="Switzerland - Males",ylim=ylimdx,
     xlab="",ylab="",axes=F,cex.main=cex.main)
axis(1)
axis(2,at=seq(0,6000,1000),labels = rep("",7))
grid();box()
points(ages,lifetable.mx(x=ages,mx=exp(LMXobs_CHEmales[,n]),sex="F")$dx,pch=1,col=col.obs[1])
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_STAD_CHEmales[,n.fore]),sex="F")$dx,col=col.mod[1],lwd=2)
lines(ages.lc_CHEmales,lifetable.mx(x=ages.lc_CHEmales,mx=exp(LMX_LC_CHEmales[,n.fore]),sex="F")$dx,col=col.mod[2],lwd=2,lty=2)
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_CODA_CHEmales[1:m,n.fore]),sex="F")$dx,col=col.mod[3],lwd=2,lty=4)
lines(ages,lifetable.mx(x=ages,mx=exp(LMX_HU_CHEmales[,n.fore]),sex="F")$dx,col=col.mod[4],lwd=2,lty=5)
legend("topleft",c("Obs, 2016","3C-STAD","LC","CODA","HU"),col=c(col.obs[1],col.mod),lty=c(NA,1,2,4,5),bty="n",
       cex=cex.leg,lwd=2,pch=c(1,NA,NA,NA,NA))
par(mfrow=c(1,1))
title(xlab = "Ages",cex.lab=cex.x.lab,
      outer = TRUE, line = 1)
dev.off()


