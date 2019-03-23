#### 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL ####
#### 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL ####
##
##  This R file contains the main function to fit the SSE 
##  model of Camarda et al. (2016)
##  
##  Authors: Ugofilippo Basellini
##           Giancarlo Camarda
##  Last update: 20/03/2019
##
## --------------- #### ------------- #### ------------- ####

# LMX_SSE=SSEfit$LHAZhat[,years%in%years.fit]
# LHAZ1=SSEfit$LHAZ1hat[,years%in%years.fit]
# LHAZ2=SSEfit$LHAZ2hat[,years%in%years.fit]
# LHAZ3=SSEfit$LHAZ3hat[,years%in%years.fit]
# years=years.fit
# E=E.fit
# Z=Z.fit
# xCHI=SSEfit$xCHI
# xSEN=SSEfit$xSEN
# xADU=SSEfit$xADU
# fore.year=years.fore[n.fore]
# stad.boot=stad.boot
# PLOT=F
# BOOT=T
# printFIT=T
# var.type_mass = "both"
# var.type_sen="trend"
# var.type_adu = "none"

# LMX_SSE=SSEfit$LHAZhat
# xCHI=SSEfit$xCHI
# LHAZ1=SSEfit$LHAZ1hat
# xSEN=SSEfit$xSEN
# LHAZ2=SSEfit$LHAZ2hat
# xADU=SSEfit$xADU
# LHAZ3=SSEfit$LHAZ3hat
# fore.year=years.fore[n.fore]
# stad.boot=stad.boot
# PLOT=F
# BOOT=T
# CENTRAL=T
# printFIT=T

## function to fit and forecast the 3C-STAD model
FitFore_STAD3C <- function(ages,years,E,Z,sex="Males",
                           LMX_SSE,xCHI,xSEN,xADU,
                           LHAZ1,LHAZ2,LHAZ3,
                           fore.year,stad.boot=10,
                           CENTRAL=T,BOOT=F,
                           PLOT=T,printFIT=F){
  
  ## loading packages
  require(MortalitySmooth)
  require(compositions)
  require(vars)
  
  ## actual hazard
  LHAZact <- log(Z/E)
  
  ## dimensions
  m <- length(ages)
  n <- length(years)
  x_chi <- xCHI
  m_chi <- length(x_chi)
  x_sen <- xSEN
  m_sen <- length(x_sen)
  x_adu <- xADU
  m_adu <- length(x_adu)
  delta <- 0.1
  delta1 <- 1
  deg <- 3
  n.simul <- stad.boot    ## increase for smoother PI
  
  ## forecast horizon
  years.fore <- (years[n]+1):fore.year
  n.fore <- length(years.fore)
  
  ## decompose FX
  FX_SSE <- FX1_SSE <- FX2_SSE <- FX3_SSE <- matrix(NA,m,n)
  i <- 1
  for (i in 1:n){
    fitFX <- FX_components(x=ages,mx = exp(LMX_SSE[,i]),x1=x_chi,x2=x_sen,x3=x_adu,
                           mx1 = exp(LHAZ1[,i]),mx2=exp(LHAZ2[,i]),mx3=exp(LHAZ3[,i]),
                           sex=sex)
    FX_SSE[,i] <- fitFX$fx
    FX1_SSE[,i] <- fitFX$fx1
    FX2_SSE[,i] <- fitFX$fx2
    FX3_SSE[,i] <- fitFX$fx3
  }
  
  ## -- CHILDHOOD ---------- 
  if (printFIT) cat("fitting Childhood component", "\n")
  
  ## actual exposures and fitted SSE deaths
  E1 <- E[ages%in%x_chi,] 
  Z1.hat <- E1*exp(LHAZ1)

  ## log death rates and densities
  LMX_chi <- LHAZ1
  FX1_chi <- FX1_SSE[ages%in%x_chi,]

  ## aligned distributions (unchanged in childhood mortality)
  LFXalign <- log(FX1_chi)
  
  ## 2. Standard = mean of the aligned density
  FXlogallmean <- exp(apply(LFXalign, 1, mean,na.rm=T))
  FXstand_chi <- FXlogallmean
  
  ## Compute hazard of standard (1y)
  ndx_chi <- floor(m_chi/3)
  lmxSTAND_chi <- hxDecr_from_fx(age = x_chi,fx = FXstand_chi, 
                                 ndx = ndx_chi,opt.print = F)$eta.hat
  
  ## 3. Coefficients of the standard 
  ## augment axis
  ages.add.r_chi <- 225
  ## define new age axis (standard needs to be expanded)
  xA_chi <- c(x_chi, seq(from=x_chi[m_chi]+delta1, by=delta1, length=ages.add.r_chi/delta1))
  mA_chi <- length(xA_chi)
  ## exclude age 0
  xA1_chi <- xA_chi[-1]
  mA1_chi <- length(xA1_chi)
  ## new B-splines parameters for standard
  xlA1_chi <- min(xA1_chi)
  xrA1_chi <- max(xA1_chi)
  ndxA1_chi <- floor((xA1_chi[mA1_chi]-xA_chi[1])/3)
  
  ## B-splines
  BA1_chi <- MortSmooth_bbase(xA1_chi, xlA1_chi, xrA1_chi, ndx=ndxA1_chi, deg=deg)
  BA_chi <- matrix(0,nrow = nrow(BA1_chi)+1, ncol = ncol(BA1_chi)+1)
  BA_chi[1,1] <- 1
  BA_chi[-1,-1] <- BA1_chi
  nbxA_chi <- ncol(BA_chi)
  fA <- c(FXstand_chi, rep(999, ages.add.r_chi/delta1)) 
  w <- c(rep(1,m_chi),rep(0,ages.add.r_chi/delta1))
  D <- diff(diag(nbxA_chi), diff=2)
  tDD <- t(D)%*%D
  P <- 10^-2 * tDD
  
  ## finding betas for the standard distribution expanded
  betasA_chi <- as.vector(solve(t(BA_chi)%*%(w*BA_chi)+P, t(BA_chi)%*%(w*log(fA))))
  coeffStand_chi <- betasA_chi
  
  ## MLE estimation ## MLE estimation
  
  ## estimate bU
  bsUP_chi <- numeric(n)
  LMXstad_chi <- matrix(NA, m_chi, n)
  CHImass <- numeric(n)
  
  i <- 1
  ## MLE estimation
  for(i in 1:n){
    ## optimization
    CHImass[i] <- sum(FX1_chi[,i])
    opt <- optimize(f=CHILD_obj_FUN,interval = c(0.1,2),x=x_chi,xA=xA_chi,
                    xup=x_chi,eta.start=lmxSTAND_chi,Dx=Z1.hat[,i], Ex=E1[,i],
                    coeff.stand=coeffStand_chi,sum.fx=CHImass[i],
                    xmin=xlA1_chi, xmax=xrA1_chi,
                    ndx=ndxA1_chi, deg=deg)
    
    ## assign 
    bUhat <- opt$minimum
    wU <- bUhat*xA_chi
    ## unique transformation function
    wb <- c(wU)
    wb1 <- wb[-1]
    ## B-splines on transformed ages
    Bwb1 <- MortSmooth_bbase(x=c(wb1),
                             xlA1_chi,xrA1_chi,ndx=ndxA1_chi,deg=deg)
    Bwb <- matrix(0,nrow = nrow(Bwb1)+1, ncol = ncol(Bwb1)+1)
    Bwb[1,1] <- 1
    Bwb[-1,-1] <- Bwb1
    ## transformed density
    fwb <- as.vector(exp(Bwb%*%coeffStand_chi))
    fwb <- fwb[xA_chi%in%x_chi]
    fwb <- fwb*CHImass[i]/sum(fwb)
    ## hazard
    eta <- hxDecr_from_fx(age = x_chi,fx=fwb,ndx=ndxA1_chi,eta.start = lmxSTAND_chi)$eta.hat
    ## assign parameter, density and compute mx
    bsUP_chi[i] <- bUhat
    LMXstad_chi[,i] <- eta
    ## plotting
    if (PLOT==TRUE){
      plot(x_chi,LMX_chi[,i],pch=16,main=paste("Childhood -",years[i]),
           ylim = range(-20,LMX_chi[1,],finite=T),xlim=range(x_chi))
      lines(x_chi,LMXstad_chi[,i],col=4,lwd=2,lty=1)
      legend("topright",c("SSE","STAD"),pch = c(16,NA),
             lty = c(0,1),col = c(1,4),bty="n",cex = 1.3,lwd = 3)
    }
  }
  
  ## -- SENESCENT ----------
  if (printFIT) cat("fitting Senescent component", "\n")
  ## original and expanded ages
  xo <- x_sen
  mo <- length(xo)
  x_sen <- seq(xo[1],120,1)
  m_sen <- length(x_sen)
  ## expanded ages at finer grid ages
  xs_sen <- seq(x_sen[1],x_sen[m_sen],delta)
  ms_sen <- length(xs_sen)
  
  ## actual exposures & fitted deaths
  E2 <- E[ages%in%xo,]
  Z2.hat <- E2*exp(LHAZ2)
  LMX_sen <- LHAZ2
  FX2_sen <- FX2_SSE[ages%in%xo,]
  
  ## B-splines parameters
  xl_sen <- min(x_sen)
  xr_sen <- max(x_sen)
  xmin_sen <- round(xl_sen - 0.01 * (xr_sen - xl_sen),3)
  xmax_sen <- round(xr_sen + 0.01 * (xr_sen - xl_sen),3)
  ndx_sen <- floor(m/3)
  
  ## B-splines bases
  B <- MortSmooth_bbase(x_sen, xmin_sen, xmax_sen, ndx_sen, deg=3)
  Bs <- MortSmooth_bbase(xs_sen, xmin_sen, xmax_sen, ndx_sen, deg=3)
  
  ## smooth mortality
  LMXs_sen <- matrix(NA, ms_sen, n)
  W_sen <- matrix(1,mo,n)
  W_sen[E2==0] <- 0
  W_senA <- rbind(W_sen,matrix(0,m_sen-mo,n))
  i <- 1
  for(i in 1:n){
    z2A <- c(Z2.hat[,i],rep(99,10))
    e2A <- c(E2[,i],rep(99,10))
    lmx.smooth.coeff <- Mort1Dsmooth(x=x_sen,y = z2A,offset = log(e2A),
                                     w = W_senA[,i],ndx = ndx_sen,deg = deg)
    LMXs_sen[,i] <- Bs %*% lmx.smooth.coeff$coefficients
  }
  
  ## compute density from smooth mortality rates
  FXs_sen <- matrix(0, nrow=ms_sen, ncol=n)
  for(i in 1:n){
    FXs_sen[,i] <- dx_from_mx(age=xs_sen,mx=exp(LMXs_sen[,i]))
  }
  
  # Mode
  M_sen <- xs_sen[apply(FXs_sen, 2, which.max)] + delta/2
  
  ## shifting parameter 
  Mstand_sen <- M_sen[1]
  s_sen <- M_sen - Mstand_sen
  
  ## aligned distributions
  FXsAlign_sen <- matrix(0, nrow=ms_sen, ncol=n)
  for(i in 1:n){
    FXsAlign_sen[,i] <- fx_shift(age=xs_sen,fx=FXs_sen[,i],shift=-s_sen[i],ndx = ndx_sen,deg = deg)
  }
  FXallmean <- exp(apply(log(FXsAlign_sen), 1, mean, na.rm=T))
  FXstand_sen <- FXallmean
  FXstand1y_sen <- FXstand_sen[xs_sen%in%x_sen]

  ## Coefficients of the standard 
  ## augment axis
  ages.add.l_sen <- 125
  ages.add.r_sen <- 150
  ## define new age axis (standard needs to be expanded)
  xA_sen <- c(rev(seq(from=x_sen[1]-delta1, by=-delta1,length=ages.add.l_sen/delta1)), 
              x_sen, 
              seq(from=x_sen[m_sen]+delta1, by=delta1, length=ages.add.r_sen/delta1))
  mA_sen <- length(xA_sen)
  
  ## new B-splines parameters for standard
  xlA_sen <- min(xA_sen)
  xrA_sen <- max(xA_sen)
  xminA_sen <- round(xlA_sen - 0.01 * (xrA_sen - xlA_sen),3)
  xmaxA_sen <- round(xrA_sen + 0.01 * (xrA_sen - xlA_sen),3)
  ndxA_sen <- floor((xA_sen[mA_sen]-xA_sen[1])/3)
  
  ## B-splines
  BA_sen <- MortSmooth_bbase(xA_sen, xminA_sen, xmaxA_sen, ndx=ndxA_sen, deg=deg)
  nbxA_sen <- ncol(BA_sen)
  
  ## coeff standards
  Standard <- coeff_stand(age=x_sen,fx=FXstand1y_sen,ndx=ndxA_sen,deg=deg,
                          ages.add.l=ages.add.l_sen,ages.add.r=ages.add.r_sen)
  coeffStand_sen <- Standard$betasA

  ## Find the break point of the age axis for each year 
  PLO_sen <- PUP_sen <-  list()
  XLO_sen <- XUP_sen <- list()
  for(i in 1:n){
    PLO_sen[[i]] <- which(xA_sen<=floor(M_sen[i]))
    PUP_sen[[i]] <- which(xA_sen>floor(M_sen[i]))
    XLO_sen[[i]] <- xA_sen[PLO_sen[[i]]]
    XUP_sen[[i]] <- xA_sen[PUP_sen[[i]]]
  }
  
  ## estimate bL and bU via MLE
  ss_sen <- bsLO_sen <- bsUP_sen <- numeric(n)
  LMXstad_sen <- matrix(NA, m_sen, n)
  SENmass <- numeric(n)
  
  i <- 1
  ## MLE estimation
  for(i in 1:n){
    conv.stad <- FALSE
    ## starting values
    if (i==1){
      start.value <- c(1,1)
    }else{
      start.value <- c(bsLO_sen[i-1],bsUP_sen[i-1])
    }
    SENmass[i] <- sum(FX2_sen[,i])
    ## MLE
    opt <- optim(par=start.value, fn=SENESCENT_obj_FUN,x=x_sen,xA=xA_sen,
                 Mstand=Mstand_sen, shat=s_sen[i],
                 xlo=XLO_sen[[i]], xup=XUP_sen[[i]],sum.fx=SENmass[i],
                 coeff.stand=coeffStand_sen, Dx=Z2.hat[,i], Ex=E2[,i],
                 xmin=xminA_sen, xmax=xmaxA_sen,
                 ndx=ndxA_sen, deg=deg)
    
    if (opt$convergence != 0) break
    
    ## assign 
    shat <- s_sen[i]
    bLhat <- opt$par[1]
    bUhat <- opt$par[2]
    
    ## compute dx
    ## segment a linear transformation function
    ## below the mode
    aL <- Mstand_sen - bLhat * (shat+Mstand_sen)
    wL <- aL + bLhat*XLO_sen[[i]]
    ## above the mode
    aU <- Mstand_sen - bUhat * (shat+Mstand_sen)
    wU <- aU + bUhat*XUP_sen[[i]]
    ## unique transformation function
    wb <- c(wL, wU)
    ## B-splines on transformed ages
    Bwb <- MortSmooth_bbase(x=c(wb),xminA_sen,xmaxA_sen,
                            ndx=ndxA_sen,deg=deg)
    ## transformed density
    fwb <- as.vector(exp(Bwb%*%coeffStand_sen))
    fwb <- fwb[xA_sen%in%x_sen]
    fwb <- fwb*SENmass[i]/sum(fwb)
    ## hazard
    mu <- convertFx(x=x_sen,data=fwb,from = "dx",to="mx")
    eta <- log(mu)
    
    ## assign parameter, density and compute mx
    ss_sen[i] <- shat
    bsLO_sen[i] <- bLhat
    bsUP_sen[i] <- bUhat
    LMXstad_sen[,i] <- eta
    ## plotting
    if (PLOT==TRUE){
      plot(xo,LMX_sen[,i],lwd=2,main=paste("Senescent -", years[i]),
           ylim = range(LMX_sen,1,finite=T),
           xlim=range(x_sen),t="p",pch=16)
      lines(x_sen,LMXstad_sen[,i],col=4,lwd=2,lty=1)
      legend("topleft",c("SSE","STAD"),pch = c(16,NA),
             lty = c(NA,1),col = c(1,4),bty="n",cex = 1.3,lwd = 3)
    }
  }
  
  ## -- EARLY-ADULTHOOD ----------
  if (printFIT) cat("fitting Early-adulthood component", "\n")
  
  ## expanded ages at finer grid ages
  xs_adu <- seq(x_adu[1],x_adu[m_adu],delta)
  ms_adu <- length(xs_adu)
  ndx_adu <- floor(m_adu/3)
  
  ## actual exposures & fitted deaths
  E3 <- E[ages%in%x_adu,]
  Z3.hat <- E3*exp(LHAZ3)
  LMX_adu <- LHAZ3
  FX3_adu <- FX3_SSE[ages%in%x_adu,]
  
  ## Standard: use the simple mean
  FXallmean <- exp(apply(log(FX3_adu), 1, mean,na.rm=T))
  FXstand1y_adu <- FXallmean

  ## Modes and shifting parameters
  M_adu <- x_adu[apply(FX3_adu, 2, which.max)] + delta1/2
  Mstand_adu <- x_adu[which.max(FXstand1y_adu)] + delta1/2
  s_adu <- M_adu - Mstand_adu

  ## Compute HAZARD of standard
  lmxSTAND_adu <- hxDecr_from_fx(age = x_adu,fx = FXstand1y_adu, ndx = ndx_adu,
                                 opt.print = FALSE)$eta.hat

  ## Coefficients of the standard 
  ## augment axis
  ages.add.l_adu <- 75
  ages.add.r_adu <- 125
  ## define new age axis (standard needs to be expanded)
  xA_adu <- c(rev(seq(from=x_adu[1]-delta1, by=-delta1,length=ages.add.l_adu/delta1)), 
              x_adu, 
              seq(from=x_adu[m_adu]+delta1, by=delta1, length=ages.add.r_adu/delta1))
  mA_adu <- length(xA_adu)
  
  ## new B-splines parameters for standard
  xlA_adu <- min(xA_adu)
  xrA_adu <- max(xA_adu)
  xminA_adu <- round(xlA_adu - 0.01 * (xrA_adu - xlA_adu),3)
  xmaxA_adu <- round(xrA_adu + 0.01 * (xrA_adu - xlA_adu),3)
  ndxA_adu <- floor((xA_adu[mA_adu]-xA_adu[1])/3)
  
  ## B-splines
  BA_adu <- MortSmooth_bbase(xA_adu, xminA_adu, xmaxA_adu, ndx=ndxA_adu, deg=deg)
  nbxA_adu <- ncol(BA_adu)
  
  ## coeff standards
  Standard <- coeff_stand(age=x_adu,fx=FXstand1y_adu,ndx=ndxA_adu,deg=deg,
                          ages.add.l=ages.add.l_adu,ages.add.r=ages.add.r_adu)
  coeffStand_adu <- Standard$betasA
  
  ## estimate b via MLE
  ss_adu <- bs_adu <- numeric(n)
  LMXstad_adu <- matrix(NA, m_adu, n)
  ADUmass <- numeric(n)
  FWB_adu <- matrix(NA, m_adu, n)
  
  i <- 1
  ## MLE estimation
  for(i in 1:n){
    conv.stad <- FALSE
    
    ## ONE PARAMETER (B)
    ADUmass[i] <- sum(FX3_adu[,i])
    ## optimization
    opt <- optimize(f=ADULT_obj_FUN,interval = c(0.25,2),x=x_adu,xA=xA_adu,
                    eta.start=lmxSTAND_adu,Dx=Z3.hat[,i], shat=s_adu[i],
                    Ex=E3[,i],Mstand=Mstand_adu,
                    coeff.stand=coeffStand_adu,sum.fx=ADUmass[i],
                    xmin=xminA_adu, xmax=xmaxA_adu,
                    ndx=ndxA_adu, deg=deg)
    
    ## assign 
    shat <- s_adu[i]
    bhat <- opt$minimum
    
    ## linear transformation function (AFT model)
    a <- Mstand_adu - bhat * (shat+Mstand_adu)
    wb <- a + bhat*xA_adu
    
    ## B-splines on transformed ages
    Bwb <- MortSmooth_bbase(x=c(wb),
                            xminA_adu,xmaxA_adu,ndx=ndxA_adu,deg=deg)
    ## transformed density
    fwb <- as.vector(exp(Bwb%*%coeffStand_adu))
    fwb <- fwb[xA_adu%in%x_adu]
    fwb <- fwb*ADUmass[i]/sum(fwb)
    ## hazard
    eta <- hxDecr_from_fx(age = x_adu,fx=fwb,ndx=ndxA_adu,
                          eta.start = lmxSTAND_adu)$eta.hat
    
    ## assign parameter, density and compute mx
    ss_adu[i] <- shat
    bs_adu[i] <- bhat
    LMXstad_adu[,i] <- eta
    FWB_adu[,i] <- fwb
    
    ## plotting
    if (PLOT==TRUE){
      plot(x_adu,LMX_adu[,i],pch=16,main=paste("Early-adulthood",years[i]),
           ylim = range(0,-30,finite=T),xlim=range(x_adu))
      lines(x_adu,LMXstad_adu[,i],col=4,lwd=2,lty=1)
      legend("topright",c("SSE","STAD"),pch = c(16,NA),
             lty = c(0,1),col = c(1,4),bty="n",cex = 1.3,lwd = 3)
    }
  }
   
  ## -- FIT COMBINE ----------
  
  ## combine mortality components
  LMX_STAD <- matrix(NA,nrow=m,ncol=n)
  for (i in 1:n){
    LMX_STAD[,i] <- log(PCLM_eta(x=ages,x1=x_chi,x2=xo,x3=x_adu,mx1=exp(LMXstad_chi[,i]),
                                 mx2=exp(LMXstad_sen[x_sen%in%xo,i]),mx3=exp(LMXstad_adu[,i])))
  }

  ## compute summary measures
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Functions")
  source("LifeTableFUN.R")
  e0.stad <- g0.stad <- numeric(n)
  for (i in 1:n){
    e0.stad[i] <- lifetable.mx(x=ages,mx=exp(LMX_STAD[,i]),sex=sex)$ex[1]
    g0.stad[i] <- GINI_func(ages=ages,mx=exp(LMX_STAD[,i]),sex=sex)
  }
  
  ## -- FORECASTING 1: MASSES ------------- 
  
  ## CoDa VAR
  df_mass <- matrix(c(SENmass,ADUmass,CHImass),ncol=3)
  df_mass_comp <- acomp(df_mass)
  df_mass_comp_alr <- alr(x=df_mass_comp)
  df_mass_comp_alr_ts <- as.ts(df_mass_comp_alr)
  var_type_mass <- "trend"
  if (cou=="SWE"&sex=="Males"&fore.year==2050) var_type_mass <- "both"
  var_m1 <- VAR(y=df_mass_comp_alr_ts,p = 1,type = var_type_mass)
  fit.var_mass <- fitted(var_m1,n.ahead =n.fore)
  res.var_mass <- residuals(var_m1,n.ahead =n.fore)
  pred.var_mass <- predict(var_m1,n.ahead =n.fore)
  df_mass_fit_alr <- data.frame(SENmass=fit.var_mass[,1],
                                ADUmass=fit.var_mass[,2])
  df_mass_fit <- alrInv(z=df_mass_fit_alr,orig = df_mass)
  df_mass_fore_alr <- data.frame(SENmass=pred.var_mass$fcst$Series.1[,1],
                                 ADUmass=pred.var_mass$fcst$Series.2[,1])
  df_mass_fore <- alrInv(z=df_mass_fore_alr,orig = df_mass)
  
  pred_sen <- df_mass_fore[,1]
  pred_adu <- df_mass_fore[,2]
  pred_chi <- df_mass_fore[,3]
  
  ## BOOTSTRAP
  ## define three bootstapping matrices
  SENmassBOOT <- CHImassBOOT <-
    ADUmassBOOT <- matrix(NA,nrow=n.fore,ncol=n.simul)
  if (BOOT){
    for(i in 1:n.simul){
      sample.res1 <- sample(c(res.var_mass[,1]),size = n,replace = T)
      sample.res2 <- sample(c(res.var_mass[,2]),size = n,replace = T)
      boot_mass <- df_mass_comp_alr_ts + cbind(sample.res1,sample.res2)
      boot_var <- VAR(y=boot_mass,p = 1,type = var_type_mass)
      boot_fore <- predict(boot_var,n.ahead = n.fore)
      boot_fore12 <- data.frame(x1=boot_fore$fcst$Series.1[,1],
                                x2=boot_fore$fcst$Series.2[,1])
      boot_fore_all <- alrInv(z=boot_fore12,orig = df_mass)
      SENmassBOOT[,i] <- boot_fore_all[,1]
      ADUmassBOOT[,i] <- boot_fore_all[,2]
      CHImassBOOT[,i] <- boot_fore_all[,3]
    }
  }
  
  ## -- FORECASTING 2: PARAMETERS ------------- 
  
  ## CHILDHOOD
  bU_chi.ts <- ts(bsUP_chi, start = years[1])
  bU_chi.mod <- auto.arima(bU_chi.ts, d=1,max.p=3,max.q=3,trace=F)
  pred.bU_chi <- forecast(bU_chi.mod, h=n.fore)
  BUchi_fore <- pred.bU_chi$mean
  
  ## bootstrap
  BUchiBOOT <- matrix(NA,nrow=n.fore,ncol=n.simul)
  if (BOOT){
    ## bootsrap S
    for(i in 1:n.simul){
      ## generate simulation with bootsrapping
      bu.sim <- simulate(bU_chi.mod, nsim=n.fore,
                         future=TRUE, bootstrap=TRUE)
      ## derive the bootsrap values
      BUchiBOOT[,i] <- bu.sim
    }
  }
  
  ## SENESCENT: S
  s_sen.ts <- ts(s_sen, start = years[1])
  s_sen.mod <- auto.arima(s_sen.ts, d=1,max.p=3,max.q=3,trace=F, ic="bic")
  pred.s_sen <- forecast(s_sen.mod, h=n.fore)
  Ssen_fore <- pred.s_sen$mean
  
  ## bootstrap
  SsenBOOT <- matrix(NA,nrow=n.fore,ncol=n.simul)
  if (BOOT){
    ## bootsrap S
    for(i in 1:n.simul){
      ## generate simulation with bootsrapping
      s.sim <- simulate(s_sen.mod, nsim=n.fore,
                        future=TRUE, bootstrap=TRUE)
      ## derive the bootsrap values
      SsenBOOT[,i] <- s.sim
    }
  }
  
  ## SENESCENT: BL, BU
  df_sen <- data.frame(bsLO_sen=bsLO_sen,bsUP_sen=bsUP_sen)
  df_sen.ts <- ts(df_sen, start = years[1])
  var_sen <- VAR(y=df_sen.ts,p = 1,type = "trend")
  var_resid_sen <- residuals(var_sen)
  var_fitted_sen <- fitted(var_sen)
  pred.var_sen <- predict(var_sen,n.ahead =n.fore)
  BLsen_fore <- pred.var_sen$fcst$bsLO_sen[,1]
  BUsen_fore <- pred.var_sen$fcst$bsUP_sen[,1]

  ## bootstrap
  BLsenBOOT <- BUsenBOOT <- matrix(NA,nrow=n.fore,ncol=n.simul)
  if (BOOT){
    for(i in 1:n.simul){
      sample.res1 <- sample(c(var_resid_sen[,1]),size = n,replace = T)
      sample.res2 <- sample(c(var_resid_sen[,2]),size = n,replace = T)
      df_sen.sim.ts <- df_sen.ts + cbind(sample.res1,sample.res2)
      var_sen_sim <- VAR(y=df_sen.sim.ts,p = 1,type = "trend")
      boot_sen_fore <- predict(var_sen_sim,n.ahead = n.fore)
      BLsenBOOT[,i] <- boot_sen_fore$fcst$bsLO_sen[,1]
      BUsenBOOT[,i] <- boot_sen_fore$fcst$bsUP_sen[,1]
    }
  }
  
  ## EARLY-ADULTHOOD
  if(all(ss_adu==0)) ss_adu <- ss_adu+rnorm(n,0,1e-1)
  var_type_adu <- "none"
  if (cou=="SWE"&sex=="Males"&fore.year==2050) var_type_adu <- "const"
  df_adu <- data.frame(ss_adu=ss_adu,bs_adu=bs_adu)
  df_adu.ts <- ts(df_adu, start = years[1])
  var_adu <- VAR(y=df_adu.ts,p = 1,type = var_type_adu)
  var_resid_adu <- residuals(var_adu)
  var_fitted_adu <- fitted(var_adu)
  pred.var_adu <- predict(var_adu,n.ahead =n.fore)
  Sadu_fore <- pred.var_adu$fcst$ss_adu[,1]
  Badu_fore <- pred.var_adu$fcst$bs_adu[,1]
  
  ## bootstrap
  SaduBOOT <- BaduBOOT <- matrix(NA,nrow=n.fore,ncol=n.simul)
  if (BOOT){
    for(i in 1:n.simul){
      sample.res1 <- sample(c(var_resid_adu[,1]),size = n,replace = T)
      sample.res2 <- sample(c(var_resid_adu[,2]),size = n,replace = T)
      df_adu.sim.ts <- df_adu.ts + cbind(sample.res1,sample.res2)
      var_adu_sim <- VAR(y=df_adu.sim.ts,p = 1,type = var_type_adu)
      boot_adu_fore <- predict(var_adu_sim,n.ahead = n.fore)
      SaduBOOT[,i] <- boot_adu_fore$fcst$ss_adu[,1]
      BaduBOOT[,i] <- boot_adu_fore$fcst$bs_adu[,1]
    }
  }
  
  ##---- COMBINE FORECASTS AND BOOTSTRAP ----
  
  ## original STAD forecast
  LMXchi_fore <- matrix(NA,m_chi,n.fore)
  LMXsen_fore <- matrix(NA,m_sen,n.fore)
  LMXadu_fore <- matrix(NA,m_adu,n.fore)
  LMX_STAD_fore <- matrix(NA,m,n.fore)
  e0stad_fore <- g0stad_fore <- rep(NA,n.fore)
  ## bootstrap results
  LMXchi_fore_boot <- matrix(NA,m_chi,n.fore)
  LMXsen_fore_boot <- matrix(NA,m_sen,n.fore)
  LMXadu_fore_boot <- matrix(NA,m_adu,n.fore)
  LMX_STAD_fore_boot <- array(NA,dim=c(m,n.fore,n.simul))
  e0stad_fore_boot <- g0stad_fore_boot <- matrix(NA,n.fore,n.simul)
  
  if (BOOT){
    cat("forecasting: computing",stad.boot,"bootstrap","\n")
    j <- 1
    k <- 1
    for (k in 1:n.simul){
      for(j in 1:n.fore) {
        ## CHILDHOOD
        ## get values of bsUP for each forecasted year
        bsUP.fore <- BUchiBOOT[j,k]
        wU <- bsUP.fore*xA_chi
        ## unique transformation function
        wb <- c(wU)
        wb1 <- wb[-1]
        ## B-splines on transformed ages
        Bwb1 <- MortSmooth_bbase(x=c(wb1),
                                 xlA1_chi,xrA1_chi,ndx=ndxA1_chi,deg=deg)
        Bwb <- matrix(0,nrow = nrow(Bwb1)+1, ncol = ncol(Bwb1)+1)
        Bwb[1,1] <- 1
        Bwb[-1,-1] <- Bwb1
        
        ## transformed density
        fwb <- as.vector(exp(Bwb%*%coeffStand_chi))
        fwb <- fwb[xA_chi%in%x_chi]
        fwb <- fwb*CHImassBOOT[j,k]/sum(fwb)
        ## hazard
        eta <- as.vector(hxDecr_from_fx(age = x_chi,
                                        fx=fwb,ndx=ndxA1_chi,eta.start = lmxSTAND_chi)$eta.hat)
        LMXchi_fore_boot[,j] <- eta
        
        ## SENESCENT
        s.fore <- SsenBOOT[j,k]
        bsLO.fore <- BLsenBOOT[j,k]
        bsUP.fore <-  BUsenBOOT[j,k]
        ## compute new mode
        Mfore_sen <- s.fore + Mstand_sen
        ## divide forecast distribution in two pieces  
        PLO_sen[[n+j]] <- which(xA_sen<floor(Mfore_sen))
        PUP_sen[[n+j]] <- which(xA_sen>=floor(Mfore_sen))
        XLO_sen[[n+j]] <- xA_sen[PLO_sen[[n+j]]]
        XUP_sen[[n+j]] <- xA_sen[PUP_sen[[n+j]]]
        
        ## segment a linear transformation function
        ## below the mode of the forecast year
        aL <- Mstand_sen - bsLO.fore * (s.fore + Mstand_sen)
        wL <- aL + bsLO.fore*XLO_sen[[n+j]]
        ## above the mode of the forecast year
        aU <- Mstand_sen - bsUP.fore * (s.fore + Mstand_sen)
        wU <- aU + bsUP.fore*XUP_sen[[n+j]]
        ## unique transformation function
        wbhat <- c(wL, wU)
        ## unique transformation function
        wb <- sort(wbhat)
        ## B-splines on transformed ages
        Bwb <- MortSmooth_bbase(x=c(wb),
                                xminA_sen,xmaxA_sen,ndx=ndxA_sen,deg=deg)
        ## transformed density
        fwb <- as.vector(exp(Bwb%*%coeffStand_sen))
        fwb <- fwb[xA_sen%in%x_sen]
        fwb <- fwb*SENmassBOOT[j,k]/sum(fwb)
        ## hazard
        mu <- convertFx(x=x_sen,data=fwb,from = "dx",to="mx")
        eta <- log(mu)
        
        ## save
        LMXsen_fore_boot[,j] <- eta
        
        ## EARLY-ADULTHOOD
        s.fore <- SaduBOOT[j,k]
        bs.fore <- BaduBOOT[j,k]
        ## linear transformation function (AFT model)
        a <- Mstand_adu - bs.fore * (s.fore+Mstand_adu)
        wb <- a + bs.fore*xA_adu
        ## B-splines on transformed ages
        Bwb <- MortSmooth_bbase(x=c(wb),
                                xminA_adu,xmaxA_adu,ndx=ndxA_adu,deg=deg)
        ## transformed density
        fwb <- as.vector(exp(Bwb%*%coeffStand_adu))
        fwb <- fwb[xA_adu%in%x_adu]
        fwb <- fwb*ADUmassBOOT[j,k]/sum(fwb)
        ## hazard
        eta <- hxDecr_from_fx(age = x_adu,fx=fwb,ndx=ndxA_adu,
                              eta.start = lmxSTAND_adu)$eta.hat
        ## save
        LMXadu_fore_boot[,j] <- eta
        
        ## put together and compute LE, gini
        LMX_STAD_fore_boot[,j,k] <- log(PCLM_eta(x=ages,x1=x_chi,x2=xo,x3=x_adu,
                                                 mx1=exp(LMXchi_fore_boot[,j]),mx2=exp(LMXsen_fore_boot[x_sen%in%xo,j]),mx3=exp(LMXadu_fore_boot[,j])))
        e0stad_fore_boot[j,k] <- lifetable.mx(x=ages,mx=exp(LMX_STAD_fore_boot[,j,k]),sex=substr(sex,1,1))$ex[1]
        g0stad_fore_boot[j,k] <- GINI_func(ages=ages,mx=exp(LMX_STAD_fore_boot[,j,k]),sex=substr(sex,1,1))
      } 
    }  
  }
  
  ## check if unplausible values
  if (any(e0stad_fore_boot<e0.stad[1],na.rm=T)){
    which.sim <- unique(which(e0stad_fore_boot<e0.stad[1],arr.ind = T)[,2])
    e0stad_fore_boot[,which.sim] <- NA
    g0stad_fore_boot[,which.sim] <- NA
    LMX_STAD_fore_boot[,,which.sim] <- NA
  }
  if (any(g0stad_fore_boot>g0.stad[1],na.rm=T)){
    which.sim <- unique(which(g0stad_fore_boot<g0.stad[1],arr.ind = T)[,2])
    e0stad_fore_boot[,which.sim] <- NA
    g0stad_fore_boot[,which.sim] <- NA
    LMX_STAD_fore_boot[,,which.sim] <- NA
  }
  
  if (CENTRAL){
    j <- 1
    for(j in 1:n.fore) {
      # cat("forecasting year",years.fore[j],"\n")
      ## CHILDHOOD
      ## get values of bsUP for each forecasted year
      bsUP.fore <- BUchi_fore[j]
      wU <- bsUP.fore*xA_chi
      ## unique transformation function
      wb <- c(wU)
      wb1 <- wb[-1]
      ## B-splines on transformed ages
      Bwb1 <- MortSmooth_bbase(x=c(wb1),
                               xlA1_chi,xrA1_chi,ndx=ndxA1_chi,deg=deg)
      Bwb <- matrix(0,nrow = nrow(Bwb1)+1, ncol = ncol(Bwb1)+1)
      Bwb[1,1] <- 1
      Bwb[-1,-1] <- Bwb1
      
      ## transformed density
      fwb <- as.vector(exp(Bwb%*%coeffStand_chi))
      fwb <- fwb[xA_chi%in%x_chi]
      fwb <- fwb*pred_chi[j]/sum(fwb)
      ## hazard
      eta <- as.vector(hxDecr_from_fx(age = x_chi,
              fx=fwb,ndx=ndxA1_chi,eta.start = lmxSTAND_chi)$eta.hat)
      LMXchi_fore[,j] <- eta
      # plot(eta)
      
      ## SENESCENT
      s.fore <- Ssen_fore[j]
      bsLO.fore <- BLsen_fore[j]
      bsUP.fore <-  BUsen_fore[j]
      ## compute new mode
      Mfore_sen <- s.fore + Mstand_sen
      ## divide forecast distribution in two pieces  
      PLO_sen[[n+j]] <- which(xA_sen<floor(Mfore_sen))
      PUP_sen[[n+j]] <- which(xA_sen>=floor(Mfore_sen))
      XLO_sen[[n+j]] <- xA_sen[PLO_sen[[n+j]]]
      XUP_sen[[n+j]] <- xA_sen[PUP_sen[[n+j]]]
      
      ## segment a linear transformation function
      ## below the mode of the forecast year
      aL <- Mstand_sen - bsLO.fore * (s.fore + Mstand_sen)
      wL <- aL + bsLO.fore*XLO_sen[[n+j]]
      ## above the mode of the forecast year
      aU <- Mstand_sen - bsUP.fore * (s.fore + Mstand_sen)
      wU <- aU + bsUP.fore*XUP_sen[[n+j]]
      ## unique transformation function
      wbhat <- c(wL, wU)
      ## unique transformation function
      wb <- sort(wbhat)
      ## B-splines on transformed ages
      Bwb <- MortSmooth_bbase(x=c(wb),
                              xminA_sen,xmaxA_sen,ndx=ndxA_sen,deg=deg)
      ## transformed density
      fwb <- as.vector(exp(Bwb%*%coeffStand_sen))
      fwb <- fwb[xA_sen%in%x_sen]
      fwb <- fwb*pred_sen[j]/sum(fwb)
      ## hazard
      mu <- convertFx(x=x_sen,data=fwb,from = "dx",to="mx")
      eta <- log(mu)
      ## save
      LMXsen_fore[,j] <- eta
      # plot(eta)
      
      ## EARLY-ADULTHOOD
      s.fore <- Sadu_fore[j]
      bs.fore <- Badu_fore[j]
      ## linear transformation function (AFT model)
      a <- Mstand_adu - bs.fore * (s.fore+Mstand_adu)
      wb <- a + bs.fore*xA_adu
      
      ## B-splines on transformed ages
      Bwb <- MortSmooth_bbase(x=c(wb),
                              xminA_adu,xmaxA_adu,ndx=ndxA_adu,deg=deg)
      ## transformed density
      fwb <- as.vector(exp(Bwb%*%coeffStand_adu))
      fwb <- fwb[xA_adu%in%x_adu]
      fwb <- fwb*pred_adu[j]/sum(fwb)
      ## hazard
      eta <- hxDecr_from_fx(age = x_adu,fx=fwb,ndx=ndxA_adu,
                            eta.start = lmxSTAND_adu)$eta.hat
      ## save
      LMXadu_fore[,j] <- eta
      # plot(eta)
      
      ## put together and compute LE, gini
      LMX_STAD_fore[,j] <- log(PCLM_eta(x=ages,x1=x_chi,x2=xo,x3=x_adu,
                                        mx1=exp(LMXchi_fore[,j]),mx2=exp(LMXsen_fore[x_sen%in%xo,j]),mx3=exp(LMXadu_fore[,j])))
      e0stad_fore[j] <- lifetable.mx(x=ages,mx=exp(LMX_STAD_fore[,j]),sex=substr(sex,1,1))$ex[1]
      g0stad_fore[j] <- GINI_func(ages=ages,mx=exp(LMX_STAD_fore[,j]),sex=substr(sex,1,1))
      
    } 
  }
  
  
  ## -- RETURN LIST ----------
  return.list <- list(LMX_STAD=LMX_STAD,LMX_STAD_fore=LMX_STAD_fore,LMX_STAD_fore_boot=LMX_STAD_fore_boot,
                      LMXstad_chi=LMXstad_chi,LMXstad_sen=LMXstad_sen,LMXstad_adu=LMXstad_adu,
                      ss_sen=ss_sen,Ssen_fore=Ssen_fore,
                      bsLO_sen=bsLO_sen,bs_adu=bs_adu,bsUP_chi=bsUP_chi,bsUP_sen=bsUP_sen,
                      BLsen_fore=BLsen_fore,Badu_fore=Badu_fore,Sadu_fore,
                      BUchi_fore=BUchi_fore,BUsen_fore=BUsen_fore,
                      CHImass=CHImass,SENmass=SENmass,ADUmass=ADUmass,
                      CHImassBOOT=CHImassBOOT,SENmassBOOT=SENmassBOOT,ADUmassBOOT=ADUmassBOOT,
                      e0.stad=e0.stad,g0.stad=g0.stad,e0stad_fore=e0stad_fore,g0stad_fore=g0stad_fore,
                      e0stad_fore_boot=e0stad_fore_boot,g0stad_fore_boot=g0stad_fore_boot,
                      SsenBOOT=SsenBOOT,BLsenBOOT=BLsenBOOT,
                      BUsenBOOT=BUsenBOOT,BaduBOOT=BaduBOOT,SaduBOOT)
  
  return(return.list)
}
