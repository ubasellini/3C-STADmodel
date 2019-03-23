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

## function to fit the SSE 2D model
Fit_SSE2D <- function(ages0,years,Z0,E0,sex="Males",
                      coef.st=NULL,
                      lambda.x1=10^3,lambda.x2=10^3.5,lambda.x3=10^3.5,
                      lambda.y1=10^2,lambda.y2=10^-5,lambda.y3=10^2.5,
                      deltay2=2,computeH=F,PRINT=T){
  ## load packages
  require(MortalitySmooth)
  require(spam)
  ## remove NA deaths
  Z0[is.na(Z0)] <- 0
  LHAZact0 <- log(Z0/E0)
  ## treat age 0 differently - exclude it from the SSE
  E <- E0[-1,]
  Z <- Z0[-1,]
  ages <- ages0[-1]
  m0 <- nrow(E0)
  ## axes and dimensions
  m <- nrow(E)
  n <- ncol(E)
  mn <- m*n
  x <- 1:m
  y <- 1:n
  ## regression weights in case of zero exposures
  ## -> interpolation
  WEI <- matrix(1,m,n)
  WEI[E==0] <- 0
  wei <- c(WEI)
  ## response and offset in column-vector
  e <- c(E)
  z <- c(Z)
  ## select over which ages each component
  ## should be restricted
  ## infant component
  staINF <- 1
  endINF <- 50
  ## aging/senescence component
  staAGI <- 1
  endAGI <- m
  ## accident-hump (middle-mortality) component
  staACC <- 1
  endACC <- 80
  ## select data & create weights for each age-window
  ## infant mortality
  p1 <- c(staINF, endINF)
  x1 <- x[x>=p1[1] & x<=p1[2]]
  W1 <- matrix(as.numeric(x%in%x1), m, n)
  Z1 <- matrix(Z[which(W1==1)], ncol=n)
  E1 <- matrix(E[which(W1==1)], ncol=n)
  ## over all mortality
  p2 <- c(staAGI, endAGI)
  x2 <- x[x>=p2[1] & x<=p2[2]]
  W2 <- matrix(as.numeric(x%in%x2), m, n)
  Z2 <- matrix(Z[which(W2==1)], ncol=n)
  E2 <- matrix(E[which(W2==1)], ncol=n)
  ## accident-hump
  p3 <- c(staACC, endACC)
  x3 <- x[x>=p3[1] & x<=p3[2]]
  W3 <- matrix(as.numeric(x%in%x3), m, n)
  Z3 <- matrix(Z[which(W3==1)], ncol=n)
  E3 <- matrix(E[which(W3==1)], ncol=n)
  ## constructing the bases
  
  ## BASIS for the childhood component
  ## B-splines over ages
  degx1 <- floor(length(x1)/3)
  xl1 <- min(x)
  xr1 <- max(x)
  xmax1 <- xr1 + 0.01 * (xr1 - xl1)
  xmin1 <- xl1 - 0.01 * (xr1 - xl1)
  Bx1 <- MortSmooth_bbase(x1, xmin1, xmax1, degx1, 3)
  nbx1 <- ncol(Bx1)
  ## B-splines over years
  degy1 <- floor(n/3)
  yl1 <- min(y)
  yr1 <- max(y)
  ymax1 <- yr1 + 0.01 * (yr1 - yl1)
  ymin1 <- yl1 - 0.01 * (yr1 - yl1)
  By1 <- MortSmooth_bbase(y, ymin1, ymax1, degy1, 3)
  nby1 <- ncol(By1)
  ## final basis for the childhood component
  B1 <- kronecker(By1, Bx1)
  
  
  ## BASIS for the aging component
  ## B-splines over ages
  degx2 <- floor(m/3)
  xl2 <- min(x)
  xr2 <- max(x)
  xmax2 <- xr2 + 0.01 * (xr2 - xl2)
  xmin2 <- xl2 - 0.01 * (xr2 - xl2)
  Bx2 <- MortSmooth_bbase(x, xmin2, xmax2, degx2, 3)
  nbx2 <- ncol(Bx2)
  ## B-splines over years
  degy2 <- floor(n/deltay2)
  yl2 <- min(y)
  yr2 <- max(y)
  ymax2 <- yr2 + 0.01 * (yr2 - yl2)
  ymin2 <- yl2 - 0.01 * (yr2 - yl2)
  By2 <- MortSmooth_bbase(y, ymin2, ymax2, degy2, 3)
  nby2 <- ncol(By2)
  ## final basis for the aging component
  B2 <- kronecker(By2, Bx2)
  
  ## BASIS for the accident-hump component
  ## B-splines over ages
  degx3 <- floor(length(x3)/3)
  xl3 <- min(x3)
  xr3 <- max(x3)
  xmax3 <- xr3 + 0.01 * (xr3 - xl3)
  xmin3 <- xl3 - 0.01 * (xr3 - xl3)
  Bx3 <- MortSmooth_bbase(x3, xmin3, xmax3, degx3, 3)
  nbx3 <- ncol(Bx3)
  ## B-splines over years
  degy3 <- floor(n/3)
  yl3 <- min(y)
  yr3 <- max(y)
  ymax3 <- yr3 + 0.01 * (yr3 - yl3)
  ymin3 <- yl3 - 0.01 * (yr3 - yl3)
  By3 <- MortSmooth_bbase(y, ymin3, ymax3, degy3, 3)
  nby3 <- ncol(By3)
  ## final basis for the accident-hump component
  B3 <- kronecker(By3, Bx3)
  
  ## complete model matrix as a list
  ## in their order and with dimensions
  XX <- list(X1=B1, X2=B2, X3=B3)
  nx <- length(XX)
  nc <- unlist(lapply(XX, ncol))
  ## indicators for the coefficients of each basis
  ind2 <- cumsum(nc)
  ind1 <- c(1, ind2[1:(nx-1)]+1)
  ind <- NULL
  for(i in 1:nx){
    ind[[i]] <- ind1[i]:ind2[i]
  }
  ## indicator for the fitted values
  ## rows
  wr1 <- rep(W1[,1], n)
  wr2 <- rep(W2[,1], n)
  wr3 <- rep(W3[,1], n)
  indR <- list(wr1=wr1, wr2=wr2, wr3=wr3)
  ## columns
  wc1 <- rep(1, nc[1])
  wc2 <- rep(1, nc[2])
  wc3 <- rep(1, nc[3])
  indC <- list(wc1=wc1, wc2=wc2, wc3=wc3)
  
  ## penalties for each component
  
  ## penalty terms for the child component
  ## including monotonity over age
  Dmon1 <- kronecker(diag(nby1),
                     diff(diag(nbx1), diff=1))
  wmon1 <- rep(0, (nbx1-1)*nby1)
  Wmon1 <- diag(wmon1)
  ## smooth penalty for the age
  Dx1 <- diff(diag(nbx1), diff=2)
  tDDx1 <- t(Dx1) %*% Dx1
  Px1 <- kronecker(diag(nby1), tDDx1)
  ## smooth penalty for the year
  Dy1 <- diff(diag(nby1), diff=2)
  tDDy1 <- t(Dy1) %*% Dy1
  Py1 <- kronecker(tDDy1, diag(nbx1))
  
  
  ## penalty terms for the aging component
  ## including monotonity over age
  Dmon2 <- kronecker(diag(nby2), diff(diag(nbx2),
                                      diff=1))
  wmon2 <- rep(0, (nbx2-1)*nby2)
  Wmon2 <- diag(wmon2)
  ## smooth penalty for the age
  Dx2 <- diff(diag(nbx2), diff=2)
  tDDx2 <- t(Dx2) %*% Dx2
  Px2 <- kronecker(diag(nby2), tDDx2)
  ## smooth penalty for the year
  Dy2 <- diff(diag(nby2), diff=2)
  tDDy2 <- t(Dy2) %*% Dy2
  Py2 <- kronecker(tDDy2, diag(nbx2))
  
  ## penalty stuff for the accident-hump component
  ## including log-concaveness
  Dcon3 <- kronecker(diag(nby3), diff(diag(nbx3),
                                      diff=2))
  wcon3 <- rep(0, (nbx3-2)*nby3)
  Wcon3 <- diag(wcon3)
  ## smooth penalty for the age
  Dx3 <- diff(diag(nbx3), diff=3)
  tDDx3 <- t(Dx3) %*% Dx3
  Px3 <- kronecker(diag(nby3), tDDx3)
  ## smooth penalty for the years
  Dy3 <- diff(diag(nby3), diff=2)
  tDDy3 <- t(Dy3) %*% Dy3
  Py3 <- kronecker(tDDy3, diag(nbx3))

  
  ## smoothing each component
  Pxy1 <- lambda.x1 * Px1 + lambda.y1 * Py1
  Pxy2 <- lambda.x2 * Px2 + lambda.y2 * Py2
  Pxy3 <- lambda.x3 * Px3 + lambda.y3 * Py3
  ## final penalty for smoothness
  P <- bdiag.spam(Pxy1, Pxy2, Pxy3)
  
  
  if (is.null(coef.st)){
    ## starting values
    
    ## fitting 2D P-splines for the senescence mortality
    ## by forcing a log-linear pattern over ages
    if (sex=="Males"){
      fit2 <- Mort2Dsmooth(x=x, y=y, Z=Z, offset=log(E),
                           W=W2*WEI,
                           method=3, lambdas=c(10^8, 1),
                           ndx = c(degx2, degy2),
                           deg = c(3,3), pord = c(2,2))
    }else{
      fit2 <- Mort2Dsmooth(x=x, y=y, Z=Z, offset=log(E),
                           W=W2*WEI,
                           method=3, lambdas=c(10^6, 1),
                           ndx = c(degx2, degy2),
                           deg = c(3,3), pord = c(2,2))
      
    }
    ## starting coefficients for the senescence component
    coef2.st <- c(fit2$coef)
    ## starting fitted values for the senescence component
    z2.st <- exp(XX[[2]] %*% coef2.st)*e ## in vector
    Z2.st <- matrix(z2.st, m, n) ## in matrix
    
    ## fitting 2D P-splines for the childhood ages
    ## taking only the first 8 ages and extrapolating
    ## for the successive ages
    WW1 <- matrix(0,length(x1),n)
    WW1[x1<=8,] <- 1
    fit1 <- Mort2Dsmooth(x=x1, y=y, Z=Z1,
                         offset=log(E1),
                         W=WW1,
                         method=3,
                         lambdas=c(10^8, 10^6),
                         ndx = c(degx1, degy1),
                         deg = c(3,3), pord = c(2,2))
    ## starting coefficients for the childhood component
    coef1.st <- c(fit1$coef)
    ## starting fitted values for the childhood component
    z1.st <- exp(XX[[1]] %*% coef1.st)*c(E1) ## in vector
    Z1.st <- matrix(z1.st, length(x1), n) ## in matrix
    ## in matrix over all ages
    Z1.stA <- matrix(0, m, n) 
    Z1.stA[W1==1] <- z1.st
    
    
    ## fitting P-splines for the accident component
    ## by fitting the differences between the
    ## overall mortality and the infant+aging components
    ## and over ages 15:40 and extrapolating
    ## previous and successive ages
    Z3.st <- Z - Z1.stA - Z2.st
    Z3.st <- matrix(Z3.st[W3==1], ncol=n)
    Z3.st[Z3.st<0] <- 0
    WW3 <- Z3.st*0
    WW3[15:40, ] <- 1
    fit3 <- Mort2Dsmooth(x=x3, y=y, Z=Z3.st,
                         offset=log(E3),
                         W=WW3,
                         method=3,
                         lambdas=c(10^4, 10^6),
                         ndx = c(degx3, degy3),
                         deg = c(3,3), pord = c(3,2))
    
    ## starting coefficients for the childhood component
    coef3.st <- c(fit3$coef)
    ## starting fitted values for the accident-hump component
    z3.st <- exp(XX[[3]] %*% coef3.st)*c(E3) ## in vector
    Z3.st <- matrix(z3.st, length(x3), n) ## in matrix
    ## in matrix over all ages
    Z3.stA <- matrix(0, m, n) 
    Z3.stA[W3==1] <- z3.st
    
    
    ## starting values at the log-mortality scale
    LHAZ1.st <- log(Z1.st/E1)
    LHAZ2.st <- log(Z2.st/E2)
    LHAZ3.st <- log(Z3.st/E3)
    Z.st <- Z1.stA + Z2.st + Z3.stA
    LHAZ.st <- log(Z.st/E)
    
    ## concatenating starting coefficients
    coef.st <- as.vector(c(coef1.st,
                           coef2.st,
                           coef3.st))
  }
  
  ## set starting coefficients
  coef <- coef.st
  
  ## shape parameter
  kappa <- 10^6
  ## solving step
  d <- 0.5
  ## small ridge penalty for numerical stability
  Pr <- 10^-4 * diag.spam(length(coef.st))
  ## maximum number of iteration
  max.it <- 200
  
  ## iteration
  for(it in 1:max.it){
    ## penalty for the monotonicity for child
    Pmon10 <- t(Dmon1) %*% Wmon1 %*% Dmon1
    Pmon1 <- kappa * Pmon10
    ## penalty for the monotonicity for aging
    Pmon20 <- t(Dmon2) %*% Wmon2 %*% Dmon2
    Pmon2 <- kappa * Pmon20
    ## penalty for the concaveness for hump
    Pcon30 <- t(Dcon3) %*% Wcon3 %*% Dcon3
    Pcon3 <- kappa * Pcon30
    ## final shape penalty term
    Psha <- bdiag.spam(Pmon1, Pmon2, Pcon3)
    ## linear predictor
    eta <- numeric(nx*mn)
    for(i in 1:nx){
      eta0 <- rep(0, mn)
      eta0[which(indR[[i]]==1)] <- XX[[i]]%*%coef[ind[[i]]]
      eta[1:mn+(i-1)*mn] <- eta0
    }
    ## components
    gamma <- exp(eta) * c(unlist(indR))
    ## expected values
    mu <- numeric(mn)
    for(i in 1:nx){
      mu <- (e * gamma[1:mn+(i-1)*mn]) + mu
    }    
    ## weights for the IWLS
    w <- mu
    ## modified model matrix for a CLM
    U <- matrix(NA, mn, sum(nc))
    for(i in 1:nx){
      u <- gamma[1:mn+(i-1)*mn]/mu * e
      u[wei==0] <- 0
      XXi <- matrix(0, nrow=mn, ncol=nc[i])
      XXi[which(indR[[i]]==1),which(indC[[i]]==1)]<-XX[[i]]
      U0 <- u * XXi
      U[,ind[[i]]] <- U0
    }
    ## regression parts for the P-CLM
    tUWU <- t(U) %*% (w * U)
    tUWUpP <- tUWU + P + Psha + Pr
    r <- z - mu
    tUr <- t(U) %*% r
    ## updating coefficients with a d-step 
    coef.old <- coef
    coef <- solve(tUWUpP, tUr + tUWU%*%coef)
    coef <- d*coef.old + (1-d)*coef
    ## update weights for shape constraints
    ## infant, monotonicity
    Wmon1.old <- Wmon1
    wmon1 <- rep(1, (nbx1-1)*nby1)
    COEF <- matrix(coef[ind[[1]]], nbx1, nby1)
    diff.COEF <- apply(COEF, 2, diff)<=0
    wmon1[c(diff.COEF)] <- 0
    Wmon1 <- diag(wmon1)
    ## aging, monotonicity
    Wmon2.old <- Wmon2
    wmon2 <- rep(1, (nbx2-1)*nby2)
    COEF <- matrix(coef[ind[[2]]], nbx2, nby2)
    diff.COEF <- apply(COEF, 2, diff)>=0
    wmon2[c(diff.COEF)] <- 0
    Wmon2 <- diag(wmon2)  
    ## accident-hump, log-concaveness
    Wcon3.old <- Wcon3
    wcon3 <- rep(0, nrow(Dcon3))
    COEF <- matrix(coef[ind[[3]]], nbx3, nby3)
    diff.COEF <- apply(COEF, 2, diff, diff=2)>=0
    wcon3[c(diff.COEF)] <- 1
    Wcon3 <- diag(wcon3)
    ## convergence criterion for the coefficients
    dif.coef<-max(abs(coef.old-coef))/max(abs(coef))
    ## stopping loop at convergence
    if (PRINT) cat(it, dif.coef, "\n")
    if(dif.coef < 1e-04 & it > 4) break
  }
  
  ## compute deviance
  zz <- z
  zz[z==0] <- 10^-8
  mumu <- mu
  mumu[mu==0] <- 10^-8
  dev <- 2*sum(z * log(zz/mumu))
  ## effective dimensions
  if (computeH){ 
    H <- solve(tUWUpP, tUWU)
    diagH <- diag(H)
    ed <- sum(diagH)
    ## BIC
    bic <- dev + log(mn)*ed
  }
  ## coefficients
  coef.hat <- coef
  
  ## linear predictor for each component
  etas <- NULL
  for(i in 1:nx){
    etas[[i]] <- XX[[i]] %*% coef.hat[ind[[i]]]
  }
  eta1.hat <- etas[[1]]
  eta2.hat <- etas[[2]]
  eta3.hat <- etas[[3]]
  
  ## linear predictor in a matrix over the whole x
  ETA.hat <- matrix(NA, mn, nx)
  for(i in 1:nx){
    ETA.hat[which(indR[[i]]==1),i] <- etas[[i]]
  }
  ## linear predictor for overall mortality
  eta.hat <- log(apply(exp(ETA.hat), 1, sum, na.rm=TRUE))
  
  ## fitted values for each component
  gamma1.hat <- exp(eta1.hat)
  gamma2.hat <- exp(eta2.hat)
  gamma3.hat <- exp(eta3.hat)
  gamma.hat <- exp(eta.hat) * c(unlist(indR))
  
  ## expected values
  mu.hat <- exp(eta.hat)*e
  
  ## deviance residuals
  res1 <- sign(z - mu.hat)
  res2 <- sqrt(2 * (z * log(z/mu.hat) - z + mu.hat))
  res <- res1 * res2 ## in vector
  RES <- matrix(res, m, n) ## in matrix
  
  ## fitted log-mortality in matrices
  LHAZ.hat <- matrix(log(mu.hat/e), m, n)
  LHAZ1.hat <- matrix(eta1.hat, length(x1), n)
  LHAZ2hat <- matrix(eta2.hat, m, n)
  LHAZ3hat <- matrix(eta3.hat, length(x3), n)
  
  ## include age 0 in LHAZhat and LHAZ1hat
  LHAZhat <- rbind(c(LHAZact0[1,]),LHAZ.hat)
  LHAZ1hat <- rbind(c(LHAZact0[1,]),LHAZ1.hat)
  
  ## compute observed and actual e0, g0
  setwd("~/Documents/Demography/Work/STADall/Github/3C-STADmodel/R/Functions")
  source("LifeTableFUN.R")
  e0act <- e0hat <- g0act <- g0hat <- numeric(n)
  for (i in 1:n){
    e0act[i] <- lifetable.mx(x=ages0,mx=exp(LHAZact0[,i]),sex=sex)$ex[1]
    e0hat[i] <- lifetable.mx(x=ages0,mx=exp(LHAZhat[,i]),sex=sex)$ex[1]
    g0act[i] <- GINI_func(ages=ages0,mx=exp(LHAZact0[,i]),sex = sex)
    g0hat[i] <- GINI_func(ages=ages0,mx=exp(LHAZhat[,i]),sex = sex)
  }
  
  ## redefine age axis
  xCHI <- c(0,x1)
  RES0 <- rbind(matrix(0,1,n),RES)
  
  ## output
  out <- list(xCHI=xCHI,xSEN=x2,xADU=x3,
              LHAZhat=LHAZhat,LHAZ1hat=LHAZ1hat,LHAZ2hat=LHAZ2hat,LHAZ3hat=LHAZ3hat,
              e0act=e0act,e0hat=e0hat,g0act=g0act,g0hat=g0hat,
              coef.hat=coef.hat,RES=RES0)
  return(out)
}