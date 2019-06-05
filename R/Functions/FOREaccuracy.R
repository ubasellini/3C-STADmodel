## ------ Function to compute forecast accuracy of a model -------- ##
##  Start date:  13/02/2018
##  Last update: 26/03/2018
##  Author: Ugofilippo Basellini
##
## --------------------------------------------------------------- ##

ForecastAccuracyIndic <- function(ages,MXact,MXforeMEAN,MXforeLOW,MXforeUP,
                                  LifeExpAct,LifeExpForeMEAN,LifeExpForeLOW,LifeExpForeUP,
                                  GiniAct,GiniForeMEAN,GiniForeLOW,GiniForeUP,
                                  alpha,ndigits=3){
  ## grab dimensions
  m <- nrow(MXact)
  n <- ncol(MXact)
  mn <- m*n
  ## POINT FORECAST ACCURACY
  LifeExpMae <- round(sum(abs(LifeExpForeMEAN-LifeExpAct))/n,ndigits)
  ## compute GINI and MAE
  GiniMae <- round(sum(abs(GiniForeMEAN-GiniAct))/n,ndigits)
  ## compute LMX and MAE
  lMXact <- log(MXact)
  lMXact[MXact==0] <- NA
  lMXfore <- log(MXforeMEAN)
  lMXfore[MXact==0] <- NA
  LmxMae <- round(sum(abs(lMXact-lMXfore),na.rm = T)/mn,ndigits)
  
  ## INTERVAL FORECAST ACCURACY
  ## LE
  if (any(LifeExpForeLOW>LifeExpForeUP)){
    LEAll <- cbind(LifeExpForeLOW,LifeExpForeUP)
    LifeExpForeLOW <- apply(LEAll,1,min)
    LifeExpForeUP <- apply(LEAll,1,max)
  }
  CoverageLE <- numeric(n)
  condLE <- (LifeExpForeLOW <= LifeExpAct) & (LifeExpAct <= LifeExpForeUP)
  CoverageLE[condLE] <- 1
  EmpCovProbLE <- sum(CoverageLE)/(n)
  CovProbDevLE <- abs(EmpCovProbLE-alpha)
  ## GINI
  if (any(GiniForeLOW>GiniForeUP)){
    GiniAll <- cbind(GiniForeLOW,GiniForeUP)
    GiniForeLOW <- apply(GiniAll,1,min)
    GiniForeUP <- apply(GiniAll,1,max)
  }
  CoverageGini <- numeric(n)
  condGini <- (GiniForeLOW <= GiniAct) & (GiniAct <= GiniForeUP)
  CoverageGini[condGini] <- 1
  EmpCovProbGini <- sum(CoverageGini)/(n)
  CovProbDevGini <- abs(EmpCovProbGini-alpha)
  ## LMX
  if (any(MXforeLOW>MXforeUP)){
    for (i in 1:n){
      MXall <- cbind(MXforeLOW[,i],MXforeUP[,i])
      MXforeLOW[,i] <- apply(MXall,1,min)
      MXforeUP[,i] <- apply(MXall,1,max)
    }
  }
  CoverageLMX <- matrix(0,nrow = m,ncol=n)
  lMXforeLOW <- log(MXforeLOW)
  lMXforeUP <- log(MXforeUP)
  condLMX <- (lMXforeLOW <= lMXact) & (lMXact <= lMXforeUP)
  CoverageLMX[condLMX] <- 1
  EmpCovProbLMX <- sum(CoverageLMX)/(mn)
  CovProbDevLMX <- round(abs(EmpCovProbLMX-alpha),ndigits)
  ## output
  estimates <- list(LifeExpAct=LifeExpAct,LifeExpForeMEAN=LifeExpForeMEAN,
                    LifeExpForeUP=LifeExpForeUP,LifeExpForeLOW=LifeExpForeLOW,
                    GiniAct=GiniAct,GiniForeMEAN=GiniForeMEAN,
                    GiniForeLOW=GiniForeLOW,GiniForeUP=GiniForeUP)
  accur <- list(LifeExpMae=LifeExpMae,GiniMae=GiniMae,LmxMae=LmxMae,
                EmpCovProbLE=round(EmpCovProbLE,ndigits),
                EmpCovProbGini=round(EmpCovProbGini,ndigits),
                EmpCovProbLMX=round(EmpCovProbLMX,ndigits))
  out <- list(estimates=estimates,accuracy=accur)
  return(out)
}

## functions to compute the Dawid-Sebastiani score
DSSfun <- function(y,BOOT){
  MU_boot <- apply(BOOT,1,mean)
  SD_boot <- apply(BOOT,1,sd)
  DSS <- (y-MU_boot)^2/(SD_boot^2) + 2*SD_boot
  return(mean(DSS))
}

DSSmatFun <- function(Y,BOOT){
  MU_boot <- apply(BOOT,c(1,2),mean)
  SD_boot <- apply(BOOT,c(1,2),sd)
  DSS <- (Y-MU_boot)^2/(SD_boot^2) + 2*SD_boot
  return(mean(DSS,na.rm = T))
}