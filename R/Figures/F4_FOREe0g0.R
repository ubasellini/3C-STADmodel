##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##-- 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL --##
##
##  R code to generate Figure 4 of the chapter 
##  Forecast life expectancy / gini 
##  
##  Authors: Ugofilippo Basellini
##           Giancarlo Camarda
##  Last update: 26/03/2019
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
e0_STAD_SWEmales <- apply(E0BOOT, 1, median, na.rm=T)
g0_STAD_SWEmales <- apply(G0BOOT, 1, median, na.rm=T)
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load CODA data
name <- paste(cou,sex,"CODA.Rdata",sep="_")
load(name)
e0_CODA_SWEmales <- CoDa.ex1$ex.med
g0_CODA_SWEmales <- CoDa.ex1$gx.med
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load HU data
name <- paste(cou,sex,"HU.Rdata",sep="_")
load(name)
e0_HU_SWEmales <- apply(E0BOOT, 1, median, na.rm=T)
g0_HU_SWEmales <- apply(G0BOOT, 1, median, na.rm=T)
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load LC data (and BT data)
name <- paste(cou,sex,"LC.Rdata",sep="_")
load(name)
e0_LC_SWEmales <- apply(E0BOOT, 1, median, na.rm=T)
g0_LC_SWEmales <- apply(G0BOOT, 1, median, na.rm=T)
e0obs_SWEmales <- e0obs
g0obs_SWEmales <- g0obs
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))


##---- SWEDISH FEMALES  ------
cou <- "SWE"
sex <- "Females"

## load STAD data
setwd(DataRes)
name <- paste(cou,sex,"3C-STAD.Rdata",sep="_")
load(name)
e0_STAD_SWEfemales <- apply(E0BOOT, 1, median, na.rm=T)
g0_STAD_SWEfemales <- apply(G0BOOT, 1, median, na.rm=T)
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load CODA data
name <- paste(cou,sex,"CODA.Rdata",sep="_")
load(name)
e0_CODA_SWEfemales <- CoDa.ex1$ex.med
g0_CODA_SWEfemales <- CoDa.ex1$gx.med
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "e0_CODA_SWEfemales","g0_CODA_SWEfemales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load HU data
name <- paste(cou,sex,"HU.Rdata",sep="_")
load(name)
e0_HU_SWEfemales <- apply(E0BOOT, 1, median, na.rm=T)
g0_HU_SWEfemales <- apply(G0BOOT, 1, median, na.rm=T)
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "e0_CODA_SWEfemales","g0_CODA_SWEfemales",
              "e0_HU_SWEfemales","g0_HU_SWEfemales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load LC data (and BT data)
name <- paste(cou,sex,"LC.Rdata",sep="_")
load(name)
e0_LC_SWEfemales <- apply(E0BOOT, 1, median, na.rm=T)
g0_LC_SWEfemales <- apply(G0BOOT, 1, median, na.rm=T)
e0obs_SWEfemales <- e0obs
g0obs_SWEfemales <- g0obs
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "e0_CODA_SWEfemales","g0_CODA_SWEfemales",
              "e0_HU_SWEfemales","g0_HU_SWEfemales",
              "e0_LC_SWEfemales","g0_LC_SWEfemales",
              "e0obs_SWEfemales","g0obs_SWEfemales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))


##---- CHE MALES  ------
cou <- "CHE"
sex <- "Males"

## load STAD data
setwd(DataRes)
name <- paste(cou,sex,"3C-STAD.Rdata",sep="_")
load(name)
e0_STAD_CHEmales <- apply(E0BOOT, 1, median, na.rm=T)
g0_STAD_CHEmales <- apply(G0BOOT, 1, median, na.rm=T)
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "e0_CODA_SWEfemales","g0_CODA_SWEfemales",
              "e0_HU_SWEfemales","g0_HU_SWEfemales",
              "e0_LC_SWEfemales","g0_LC_SWEfemales",
              "e0obs_SWEfemales","g0obs_SWEfemales",
              "e0_STAD_CHEmales","g0_STAD_CHEmales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load CODA data
name <- paste(cou,sex,"CODA.Rdata",sep="_")
load(name)
e0_CODA_CHEmales <- CoDa.ex1$ex.med
g0_CODA_CHEmales <- CoDa.ex1$gx.med
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "e0_CODA_SWEfemales","g0_CODA_SWEfemales",
              "e0_HU_SWEfemales","g0_HU_SWEfemales",
              "e0_LC_SWEfemales","g0_LC_SWEfemales",
              "e0obs_SWEfemales","g0obs_SWEfemales",
              "e0_STAD_CHEmales","g0_STAD_CHEmales",
              "e0_CODA_CHEmales","g0_CODA_CHEmales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load HU data
name <- paste(cou,sex,"HU.Rdata",sep="_")
load(name)
e0_HU_CHEmales <- apply(E0BOOT, 1, median, na.rm=T)
g0_HU_CHEmales <- apply(G0BOOT, 1, median, na.rm=T)
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "e0_CODA_SWEfemales","g0_CODA_SWEfemales",
              "e0_HU_SWEfemales","g0_HU_SWEfemales",
              "e0_LC_SWEfemales","g0_LC_SWEfemales",
              "e0obs_SWEfemales","g0obs_SWEfemales",
              "e0_STAD_CHEmales","g0_STAD_CHEmales",
              "e0_CODA_CHEmales","g0_CODA_CHEmales",
              "e0_HU_CHEmales","g0_HU_CHEmales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load LC data (and BT data)
name <- paste(cou,sex,"LC.Rdata",sep="_")
load(name)
e0_LC_CHEmales <- apply(E0BOOT, 1, median, na.rm=T)
g0_LC_CHEmales <- apply(G0BOOT, 1, median, na.rm=T)
e0obs_CHEmales <- e0obs
g0obs_CHEmales <- g0obs
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "e0_CODA_SWEfemales","g0_CODA_SWEfemales",
              "e0_HU_SWEfemales","g0_HU_SWEfemales",
              "e0_LC_SWEfemales","g0_LC_SWEfemales",
              "e0obs_SWEfemales","g0obs_SWEfemales",
              "e0_STAD_CHEmales","g0_STAD_CHEmales",
              "e0_CODA_CHEmales","g0_CODA_CHEmales",
              "e0_HU_CHEmales","g0_HU_CHEmales",
              "e0_LC_CHEmales","g0_LC_CHEmales",
              "e0obs_CHEmales","g0obs_CHEmales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

##---- CHE FEMALES  ------
cou <- "CHE"
sex <- "Females"

## load STAD data
setwd(DataRes)
name <- paste(cou,sex,"3C-STAD.Rdata",sep="_")
load(name)
e0_STAD_CHEfemales <- apply(E0BOOT, 1, median, na.rm=T)
g0_STAD_CHEfemales <- apply(G0BOOT, 1, median, na.rm=T)
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "e0_CODA_SWEfemales","g0_CODA_SWEfemales",
              "e0_HU_SWEfemales","g0_HU_SWEfemales",
              "e0_LC_SWEfemales","g0_LC_SWEfemales",
              "e0obs_SWEfemales","g0obs_SWEfemales",
              "e0_STAD_CHEmales","g0_STAD_CHEmales",
              "e0_CODA_CHEmales","g0_CODA_CHEmales",
              "e0_HU_CHEmales","g0_HU_CHEmales",
              "e0_LC_CHEmales","g0_LC_CHEmales",
              "e0obs_CHEmales","g0obs_CHEmales",
              "e0_STAD_CHEfemales","g0_STAD_CHEfemales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load CODA data
name <- paste(cou,sex,"CODA.Rdata",sep="_")
load(name)
e0_CODA_CHEfemales <- CoDa.ex1$ex.med
g0_CODA_CHEfemales <- CoDa.ex1$gx.med
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "e0_CODA_SWEfemales","g0_CODA_SWEfemales",
              "e0_HU_SWEfemales","g0_HU_SWEfemales",
              "e0_LC_SWEfemales","g0_LC_SWEfemales",
              "e0obs_SWEfemales","g0obs_SWEfemales",
              "e0_STAD_CHEmales","g0_STAD_CHEmales",
              "e0_CODA_CHEmales","g0_CODA_CHEmales",
              "e0_HU_CHEmales","g0_HU_CHEmales",
              "e0_LC_CHEmales","g0_LC_CHEmales",
              "e0obs_CHEmales","g0obs_CHEmales",
              "e0_STAD_CHEfemales","g0_STAD_CHEfemales",
              "e0_CODA_CHEfemales","g0_CODA_CHEfemales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load HU data
name <- paste(cou,sex,"HU.Rdata",sep="_")
load(name)
e0_HU_CHEfemales <- apply(E0BOOT, 1, median, na.rm=T)
g0_HU_CHEfemales <- apply(G0BOOT, 1, median, na.rm=T)
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "e0_CODA_SWEfemales","g0_CODA_SWEfemales",
              "e0_HU_SWEfemales","g0_HU_SWEfemales",
              "e0_LC_SWEfemales","g0_LC_SWEfemales",
              "e0obs_SWEfemales","g0obs_SWEfemales",
              "e0_STAD_CHEmales","g0_STAD_CHEmales",
              "e0_CODA_CHEmales","g0_CODA_CHEmales",
              "e0_HU_CHEmales","g0_HU_CHEmales",
              "e0_LC_CHEmales","g0_LC_CHEmales",
              "e0obs_CHEmales","g0obs_CHEmales",
              "e0_STAD_CHEfemales","g0_STAD_CHEfemales",
              "e0_CODA_CHEfemales","g0_CODA_CHEfemales",
              "e0_HU_CHEfemales","g0_HU_CHEfemales",
              "cou","sex","DataRes")
rm(list=setdiff(ls(),var.keep))

## load LC data (and BT data)
name <- paste(cou,sex,"LC.Rdata",sep="_")
load(name)
e0_LC_CHEfemales <- apply(E0BOOT, 1, median, na.rm=T)
g0_LC_CHEfemales <- apply(G0BOOT, 1, median, na.rm=T)
e0obs_CHEfemales <- e0obs
g0obs_CHEfemales <- g0obs
var.keep <- c("e0_STAD_SWEmales","g0_STAD_SWEmales",
              "e0_CODA_SWEmales","g0_CODA_SWEmales",
              "e0_HU_SWEmales","g0_HU_SWEmales",
              "e0_LC_SWEmales","g0_LC_SWEmales",
              "e0obs_SWEmales","g0obs_SWEmales",
              "e0_STAD_SWEfemales","g0_STAD_SWEfemales",
              "e0_CODA_SWEfemales","g0_CODA_SWEfemales",
              "e0_HU_SWEfemales","g0_HU_SWEfemales",
              "e0_LC_SWEfemales","g0_LC_SWEfemales",
              "e0obs_SWEfemales","g0obs_SWEfemales",
              "e0_STAD_CHEmales","g0_STAD_CHEmales",
              "e0_CODA_CHEmales","g0_CODA_CHEmales",
              "e0_HU_CHEmales","g0_HU_CHEmales",
              "e0_LC_CHEmales","g0_LC_CHEmales",
              "e0obs_CHEmales","g0obs_CHEmales",
              "e0_STAD_CHEfemales","g0_STAD_CHEfemales",
              "e0_CODA_CHEfemales","g0_CODA_CHEfemales",
              "e0_HU_CHEfemales","g0_HU_CHEfemales",
              "e0_LC_CHEfemales","g0_LC_CHEfemales",
              "e0obs_CHEfemales","g0obs_CHEfemales",
              "cou","sex","DataRes","years","years.fore")
rm(list=setdiff(ls(),var.keep))

##---- PLOTTING -----
library("RColorBrewer")
col.obs <- c("grey20","grey40")
ylimE0_SWE <- range(e0obs_SWEmales,e0obs_SWEfemales,e0_STAD_SWEmales,e0_STAD_SWEfemales,
                    e0_CODA_SWEfemales)
ylimG0_SWE <- range(g0obs_SWEmales,g0obs_SWEfemales,g0_STAD_SWEmales,g0_STAD_SWEfemales,
                 g0_CODA_SWEfemales)
ylimE0_CHE <- range(e0obs_CHEmales,e0obs_CHEfemales,e0_STAD_CHEmales,e0_STAD_CHEfemales)
ylimG0_CHE <- range(g0obs_CHEfemales,g0obs_CHEmales,g0_STAD_CHEfemales,g0_STAD_CHEmales,
                 g0_CODA_CHEfemales)
ylimE0 <- range(ylimE0_SWE,ylimE0_CHE)
ylimG0 <- range(ylimG0_SWE,ylimG0_CHE)
col.mod <- c(brewer.pal(n=8,"Blues")[7],
             brewer.pal(n=8,"Greens")[7],
             brewer.pal(n=8,"Oranges")[6],
             brewer.pal(n=8,"Purples")[8])
cex.main <- 1.3
cex.x.lab <- 1.75
cex.y.lab <- 1.5
cex.x.axis <- 1.25
cex.y.axis <- 0.9
cex.leg <- 1.3

setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/Manuscript/Figures")
pdf("F4.pdf",width = 10, height = 8)
op <- par(mfrow = c(2,2),
          oma = c(2,2,0,0) + 0.1,
          mar = c(1,2,1,1) + 0.1)
## bottom, left, top, right
## E0 - Sweden
plot(years,e0obs_SWEfemales,t="n",main="Life Expectancy - Sweden",ylim=ylimE0,xlim=range(years,years.fore),
     xlab="",ylab="",axes=F,cex.main=cex.main)
axis(1,at=seq(1960,2040,20),labels = rep("",5))
axis(2,las=2,cex.axis=cex.y.axis)
grid();box()
mtext(text = expression(e^0),side=2,
      line = 2.25,cex=cex.y.lab,las=1)
points(years,e0obs_SWEmales,pch=1,col=col.obs[2])
points(years,e0obs_SWEfemales,pch=5,col=col.obs[1])
lines(years.fore,e0_STAD_SWEmales,col=col.mod[1],lwd=2)
lines(years.fore,e0_LC_SWEmales,col=col.mod[2],lwd=2,lty=2)
lines(years.fore,e0_CODA_SWEmales,col=col.mod[3],lwd=2,lty=4)
lines(years.fore,e0_HU_SWEmales,col=col.mod[4],lwd=2,lty=5)
lines(years.fore,e0_STAD_SWEfemales,col=col.mod[1],lwd=2)
lines(years.fore,e0_LC_SWEfemales,col=col.mod[2],lwd=2,lty=2)
lines(years.fore,e0_CODA_SWEfemales,col=col.mod[3],lwd=2,lty=4)
lines(years.fore,e0_HU_SWEfemales,col=col.mod[4],lwd=2,lty=5)
legend("topleft",c("Females","Males"),col=col.obs,pch=c(5,1),bty="n",
       cex=cex.leg)

## E0 - Switzerland
plot(years,e0obs_CHEmales,t="n",main="Life Expectancy - Switzerland",ylim=ylimE0,xlim=range(years,years.fore),
     xlab="",ylab="",axes=F,cex.main=cex.main)
axis(1,at=seq(1960,2040,20),labels = rep("",5))
axis(2,las=2,cex.axis=cex.y.axis)
grid();box()
points(years,e0obs_CHEmales,pch=1,col=col.obs[2])
points(years,e0obs_CHEfemales,pch=5,col=col.obs[1])
lines(years.fore,e0_STAD_CHEmales,col=col.mod[1],lwd=2)
lines(years.fore,e0_LC_CHEmales,col=col.mod[2],lwd=2,lty=2)
lines(years.fore,e0_CODA_CHEmales,col=col.mod[3],lwd=2,lty=4)
lines(years.fore,e0_HU_CHEmales,col=col.mod[4],lwd=2,lty=5)
lines(years.fore,e0_STAD_CHEfemales,col=col.mod[1],lwd=2)
lines(years.fore,e0_LC_CHEfemales,col=col.mod[2],lwd=2,lty=2)
lines(years.fore,e0_CODA_CHEfemales,col=col.mod[3],lwd=2,lty=4)
lines(years.fore,e0_HU_CHEfemales,col=col.mod[4],lwd=2,lty=5)
legend("bottomright",c("3C-STAD","LC","CODA","HU"),col=col.mod,lty=c(1,2,4,5),bty="n",
       cex=cex.leg,lwd=2)


## G0 - Sweden
plot(years,g0obs_SWEfemales,t="n",main="Gini - Sweden",ylim=ylimG0,xlim=range(years,years.fore),
     xlab="",ylab="",axes=F,cex.main=cex.main)
axis(1,at=seq(1960,2040,20))
axis(2,las=2,cex.axis=cex.y.axis)
grid();box()
mtext(text = expression(g^0),side=2,
      line = 2.25,cex=cex.y.lab,las=1)
points(years,g0obs_SWEmales,pch=1,col=col.obs[1])
points(years,g0obs_SWEfemales,pch=5,col=col.obs[2])
lines(years.fore,g0_STAD_SWEfemales,col=col.mod[1],lwd=2)
lines(years.fore,g0_LC_SWEfemales,col=col.mod[2],lwd=2,lty=2)
lines(years.fore,g0_CODA_SWEfemales,col=col.mod[3],lwd=2,lty=4)
lines(years.fore,g0_HU_SWEfemales,col=col.mod[4],lwd=2,lty=5)
lines(years.fore,g0_STAD_SWEmales,col=col.mod[1],lwd=2)
lines(years.fore,g0_LC_SWEmales,col=col.mod[2],lwd=2,lty=2)
lines(years.fore,g0_CODA_SWEmales,col=col.mod[3],lwd=2,lty=4)
lines(years.fore,g0_HU_SWEmales,col=col.mod[4],lwd=2,lty=5)
legend("topright",c("3C-STAD","LC","CODA","HU"),col=col.mod,lty=c(1,2,4,5),bty="n",
       cex=cex.leg,lwd=2)

## G0 - Switzerland
plot(years,g0obs_CHEmales,t="n",main="Gini - Switzerland",ylim=ylimG0,xlim=range(years,years.fore),
     xlab="",ylab="",axes=F,cex.main=cex.main)
axis(1,at=seq(1960,2040,20))
axis(2,las=2,cex.axis=cex.y.axis)
grid();box()
points(years,g0obs_CHEmales,pch=1,col=col.obs[2])
points(years,g0obs_CHEfemales,pch=5,col=col.obs[1])
lines(years.fore,g0_STAD_CHEmales,col=col.mod[1],lwd=2)
lines(years.fore,g0_LC_CHEmales,col=col.mod[2],lwd=2,lty=2)
lines(years.fore,g0_CODA_CHEmales,col=col.mod[3],lwd=2,lty=4)
lines(years.fore,g0_HU_CHEmales,col=col.mod[4],lwd=2,lty=5)
lines(years.fore,g0_STAD_CHEfemales,col=col.mod[1],lwd=2)
lines(years.fore,g0_LC_CHEfemales,col=col.mod[2],lwd=2,lty=2)
lines(years.fore,g0_CODA_CHEfemales,col=col.mod[3],lwd=2,lty=4)
lines(years.fore,g0_HU_CHEfemales,col=col.mod[4],lwd=2,lty=5)
legend("topright",c("Females","Males"),col=col.obs,pch=c(5,1),bty="n",
       cex=cex.leg)
par(mfrow=c(1,1))
title(xlab = "Year",cex.lab=cex.x.lab,
      outer = TRUE, line = 1)
dev.off()


