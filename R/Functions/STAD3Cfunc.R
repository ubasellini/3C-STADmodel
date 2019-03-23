#### 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL ####
#### 3C-STAD MODEL #### 3C-STAD MODEL #### 3C-STAD MODEL ####
##
##  This R file contains the functions called by the  to run the 
##  main function FitFore3CSTAD.R
##  
##  Authors: Ugofilippo Basellini
##           Giancarlo Camarda
##  Last update: 20/03/2019
##
## --------------- #### ------------- #### ------------- ####

##--- SENESCENT COMPONENT FUNCTIONS ------------- 

## Function to compute dx from mx 
dx_from_mx <- function(age,mx){
  ## length of age range
  xs <- age
  ms <- length(xs)
  delta <- round(diff(xs)[1],4)
  ## identity matrix and C matrix
  I <- diag(ms)
  C <- lower.tri(I, diag = TRUE)
  C[C==1] <- -delta 
  ## compute dx
  dx <- mx * exp(C%*%mx)
  return(dx)
}

## function for aligning fx given a shifting parameter
fx_shift <- function(age,fx,shift,ndx=25,deg=3){
  ## length of age range
  xs <- age
  ms <- length(xs)
  delta <- round(diff(age)[1],4)
  ages.add <- 0
  ## augment basis
  if (shift > 0) ages.add = ceiling(shift + 20)
  if (shift < 0) ages.add = ceiling(abs(shift - 20))
  ## weights
  w <- c(rep(0,ages.add/delta),rep(1,ms),rep(0,ages.add/delta))
  ## augmented data
  xA <- c(rev(seq(from=xs[1]-delta, by=-delta,length=ages.add/delta)), xs, 
          seq(from=xs[ms]+delta, by=delta, length=ages.add/delta))
  fA <- c(rep(999, ages.add/delta), fx, rep(999, ages.add/delta))
  ## new basis & B-splines parameters
  xl <- min(xA)
  xr <- max(xA)
  xmin <- round(xl - 0.01 * (xr - xl),3)
  xmax <- round(xr + 0.01 * (xr - xl),3)
  nb <- ndx+deg
  ## new B-splines bases
  BA <- MortSmooth_bbase(xA, xmin, xmax, ndx, deg)
  ## smoothing + extrapolation
  D <- diff(diag(nb), diff=2)
  tDD <- t(D)%*%D
  P <- 10^-2 * tDD
  betasA <- solve(t(BA)%*%(w*BA)+P, t(BA)%*%(w*log(fA)))
  ## shifting
  ws.i <- xs - shift
  ## evaluating B on shifted x
  Bshift <- MortSmooth_bbase(ws.i, xl=xmin, xr=xmax, ndx, deg)
  ## shifted density
  fshift <- exp(Bshift %*% betasA)
  out <- fshift
  if (shift == 0) out <- fx
  return(out)
}

## function for obtaining coefficients of the standard distribution
coeff_stand <- function(age,fx,ages.add.l=30,ages.add.r=20,ndx=30,deg=3){
  ## Extrapolate age range to left and right + derive coefficients
  ## length of age range
  xs <- age
  delta <- round(diff(xs)[1],4)
  ms <- length(xs)
  ## weights
  w <- c(rep(0,ages.add.l/delta),rep(1,ms),rep(0,ages.add.r/delta))
  ## augmented data
  xA <- c(rev(seq(from=xs[1]-delta, by=-delta,length=ages.add.l/delta)), 
          xs, 
          seq(from=xs[ms]+delta, by=delta, length=ages.add.r/delta))
  fA <- c(rep(999, ages.add.l/delta), fx, 
          rep(999, ages.add.r/delta))
  ## new basis & B-splines parameters
  xl <- min(xA)
  xr <- max(xA)
  xmin <- round(xl - 0.01 * (xr - xl),3)
  xmax <- round(xr + 0.01 * (xr - xl),3)
  nb <- ndx+deg
  ## B-splines bases standard expanded
  BA <- MortSmooth_bbase(xA, xmin, xmax, ndx, deg)
  ## simple smoothing + extrapolation
  D <- diff(diag(nb), diff=2)
  tDD <- t(D)%*%D
  P <- 10^-5 * tDD
  ## finding betas for the standard distribution expanded
  betasA <- solve(t(BA)%*%(w*BA)+P, t(BA)%*%(w*log(fA)))
  ## effective dimension
  tBWB <- t(BA)%*%(w*BA)
  tBWBpP <- tBWB + P 
  H <- solve(tBWBpP, tBWB)
  h <- diag(H)
  ed <- sum(h)
  ## return
  out <- list(betasA=betasA,ed=ed)
  return(out)
}

## optimization function for bL and bU
SENESCENT_obj_FUN <- function(par,x,xA,sum.fx,
                              Mstand, shat, xlo, xup,coeff.stand, Dx, Ex, 
                              xmin, xmax, ndx, deg){
  require(MortalityLaws)
  ## starting b
  bL <- par[1]
  bU <- par[2]
  ## segment a linear transformation function
  ## below the mode
  aL <- Mstand - bL * (shat+Mstand)
  wL <- aL + bL*xlo
  ## above the mode
  aU <- Mstand - bU * (shat+Mstand)
  wU <- aU + bU*xup
  ## unique transformation function
  wb <- c(wL, wU)
  wb <- sort(wb)
  ## B-splines on transformed ages
  Bwb <- MortSmooth_bbase(x=c(wb),
                          xmin,xmax,ndx=ndx,deg=deg)
  ## transformed density
  dwb <- as.vector(exp(Bwb%*%coeff.stand))
  dwb <- dwb[xA%in%x]
  dwb <- dwb*sum.fx/sum(dwb)
  # if (sum(dwb) > 1) fwb <- dwb/sum(dwb)
  ## hazard
  mu <- convertFx(x=x,data=dwb,from = "dx",to="mx")
  eta <- log(mu)
  # eta <- log(mx_from_dx(dx=dwb))
  # mu <- exp(eta)
  ## minimise minus the Log-Likelihood (maximise the LL)
  Lij <- -sum(Dx * eta[1:length(Dx)]-Ex*mu[1:length(Ex)], na.rm = T)
  return(Lij)
}


##--- CHILDHOOD COMPONENT FUNCTIONS ------------- 

## optimization function for Childhood mortality
CHILD_obj_FUN <- function(par,x,xA,eta.start,
                           xup,coeff.stand, Dx, Ex, sum.fx,
                           xmin, xmax, ndx, deg){
  ## starting b
  bU <- par[1]
  ## segment a linear transformation function
  ## above the mode
  wU <- bU*xup
  ## unique transformation function
  wb <- c(wU)
  wb1 <- wb[-1]
  ## B-splines on transformed ages
  Bwb1 <- MortSmooth_bbase(x=c(wb1),
                           xmin,xmax,ndx=ndx,deg=deg)
  Bwb <- matrix(0,nrow = nrow(Bwb1)+1, ncol = ncol(Bwb1)+1)
  Bwb[1,1] <- 1
  Bwb[-1,-1] <- Bwb1
  ## transformed density
  fwb <- as.vector(exp(Bwb%*%coeff.stand))
  fwb <- fwb[xA%in%x]
  if (sum(fwb) > 1) fwb <- fwb/sum(fwb)
  fwb <- fwb*sum.fx/sum(fwb)
  ## hazard
  eta <- hxDecr_from_fx(age = x,fx=fwb,ndx=ndx,eta.start = eta.start)$eta.hat
  mu <- exp(eta)
  ## minimise minus the Log-Likelihood (max LL)
  Lij <- -sum(Dx * eta-Ex*mu)
  return(Lij)
}

## Function to compute fx from a SMOOTH hx 
fx_from_hx <- function(age,hx){
  ## length of age range
  xs <- age
  ms <- length(xs)
  delta <- round(diff(xs)[1],4)
  ## identity matrix and C matrix
  I <- diag(ms)
  C <- lower.tri(I, diag = TRUE)
  C[C==1] <- -delta 
  ## compute fx
  fx <- hx * exp(C%*%hx)
  return(fx)
}

## Function to compute the gradient of fx from a SMOOTH hx 
gradfx_from_hx <- function(age,hx){
  ## length of age range
  xs <- age
  ms <- length(xs)
  delta <- round(diff(xs)[1],4)
  ## identity matrix and C matrix
  I <- diag(ms)
  C <- lower.tri(I, diag = TRUE)
  C[C==1] <- -delta 
  ## compute gradient
  v <- hx * exp(C%*%hx)
  M1a <- C %*% diag(as.vector(hx))
  M1 <- M1a + I
  V <- diag(as.vector(v))
  G <- V %*% M1
  return(G)
}

## Function to compute hx from a SMOOTH fx 
hxDecr_from_fx <- function(age,fx,ndx,eta.start=NA,max.it=40,opt.print=FALSE){
  require(MortalitySmooth)
  ## length of age range
  x <- age
  m <- length(x)
  delta <- diff(x)[1]
  ## B-splines parameters
  xl <- min(x)
  xr <- max(x)
  xmin <- round(xl - 0.01 * (xr - xl),3)
  xmax <- round(xr + 0.01 * (xr - xl),3)
  deg <- 3
  nbx <- ndx+deg
  ## B-splines bases
  B <- MortSmooth_bbase(x, xmin, xmax, ndx, deg)
  
  ## Derive ETA1 
  ## numerical approximation of survival 
  Sx <- numeric(m)
  Sx <- 1 - cumsum(fx)
  Sx[Sx==0] <- min(Sx[Sx!=0])
  ## starting numerical eta
  options(warn=-1)
  eta0 <- log(fx/Sx)
  options(warn=0)
  eta0[is.infinite(eta0)] <- 0
  eta0[is.nan(eta0)] <- 0
  
  ## Series of LAMBDAS - smoothing needed to find a solution 
  ## (start from lowest smoothing possible)
  LAMBDAS <- 10^(c(-3,0,2))
  nlambdas <- length(LAMBDAS)
  
  ## Collect starting eta in a list
  if(is.na(eta.start)[1]){
    ETAS <- list(eta0=eta0)
    # ETAS <- list(eta1=eta1,eta2=eta2,eta3=eta3)
  }else{
    ETAS <- list(eta.start=eta.start,eta0=eta0)
    # ETAS <- list(eta.start=eta.start,eta1=eta1,eta2=eta2,eta3=eta3)
  }
  netas <- length(ETAS)
  
  ## Find solution
  for (l in 1:nlambdas){
    lambda <- LAMBDAS[l]
    if (opt.print== TRUE) cat("searching lambda =", lambda, "\n")
    for (i in 1:netas){
      if (opt.print== TRUE) cat("optimizing", names(ETAS)[i], "\n")
      eta <- ETAS[[i]]
      eta.hat <- hxDecr_iter_fun(age=x,fx=fx,ndx=ndx,eta=eta,lambda=LAMBDAS[l],max.it=max.it)
      if (eta.hat$conv.flag == TRUE) break
    }
    if (eta.hat$conv.flag == TRUE) break
  }
  out <- list(eta.hat=eta.hat$eta,conv=eta.hat$conv.flag)
  return(out)
}

hxDecr_iter_fun <- function(age,fx,ndx,eta,lambda,max.it=40){
  ## dimensions
  x <- age
  m <- length(age)
  ## B-splines parameters
  xl <- min(x)
  xr <- max(x)
  xmin <- round(xl - 0.01 * (xr - xl),3)
  xmax <- round(xr + 0.01 * (xr - xl),3)
  deg <- 3
  nbx <- ndx+deg
  ## B-splines bases
  B <- MortSmooth_bbase(x, xmin, xmax, ndx, deg)
  
  ## define smoothness penalty
  D_1y <- diff(diag(m), diff=2)
  tDD_1y <- t(D_1y)%*%D_1y
  P_1y <- lambda * tDD_1y
  
  ## Convergence Flag
  conv.flag <- "FALSE"
  ## LOOP
  for(it in 1:max.it){
    ## current f
    f <- fx_from_hx(age=x,hx=exp(eta))
    ## gradient
    G <- gradfx_from_hx(age=x,hx=exp(eta))
    ## residuals
    r <- fx - f
    f2 <- as.vector(f^2)
    ## regression
    WG <- G / f2
    tGWG <- t(G) %*% WG
    tGWGpP <- tGWG + P_1y 
    Wr <- r/f2
    tGWr <- t(G) %*% Wr + tGWG %*% eta
    ## old eta
    eta.old <- eta
    ## new eta
    eta <- solve(tGWGpP, tGWr)
    ## break loop if regression goes wrong
    if (is.nan(eta[m]) | any(is.na(eta)) | eta[m] > 15) break
    ## convergence
    conv <- max(abs(eta-eta.old))
    # cat(it,conv,"\n")
    if(conv < 10^-4 & it > 4){
      if(diff(eta)[m-1]>0){
        break
      }else{
        conv.flag <- "TRUE"
        break
      }
    } 
  }
  out <- list(eta=eta,conv.flag=conv.flag)
  return(out)
}


##--- EARLY-ADULTHOOD COMPONENT FUNCTIONS ------------- 

## optimization function for Early-adulthood mortality
ADULT_obj_FUN <- function(par,x,xA,eta.start,Mstand,shat,
                             coeff.stand, Dx, Ex, sum.fx,
                             xmin, xmax, ndx, deg){
  ## starting b
  b <- par[1]
  ## linear transformation function (AFT model)
  a <- Mstand - b * (shat+Mstand)
  wb <- a + b*xA
  ## unique transformation function
  ## B-splines on transformed ages
  Bwb <- MortSmooth_bbase(x=c(wb),
                          xmin,xmax,ndx=ndx,deg=deg)
  ## transformed density
  fwb <- as.vector(exp(Bwb%*%coeff.stand))
  fwb <- fwb[xA%in%x]
  fwb <- fwb*sum.fx/sum(fwb)
  ## hazard
  eta <- hxDecr_from_fx(age = x,fx=fwb,ndx=ndx,eta.start = eta.start)$eta.hat
  mu <- exp(eta)
  ## minimise minus the Log-Likelihood (max LL)
  Lij <- -sum(Dx * eta-Ex*mu)
  return(Lij)
}

##--- OTHER GENERIC FUNCTIONS ------------- 

## function for deriving component-specific fx
FX_components <- function(x,mx,x1,mx1,x2,mx2,x3,mx3,sex="Male"){
  ## dimensions
  m <- length(x)
  m1 <- length(x1)
  m2 <- length(x2)
  m3 <- length(x3)
  ## derive LT for three hazards
  ## first way - some problems at age 110
  LT <- lifetable.mx(x=x,mx=mx,sex = substr(sex,1,1))
  fx <- (LT$dx)/(LT$lx[1])
  lx <- LT$lx
  ## exclude NaN
  fx[is.nan(fx)] <- 0
  lx[is.nan(lx)] <- 0
  ## rescale fx
  fx <- fx/sum(fx)
  # plot(x,log(fx),t="l",lwd=2)
  fx1 <- fx2 <- fx3 <- numeric(m)
  fx1[x%in%x1] <- lx[x%in%x1]*mx1
  fx2[x%in%x2] <- lx[x%in%x2]*mx2
  fx3[x%in%x3] <- lx[x%in%x3]*mx3
  ## adjust at each age to sum to fx
  i <- m
  for (i in 1:m){
    fxTot <- fx[i]
    if (fxTot!=0){
      fxsTot <- sum(c(fx1[i],fx2[i],fx3[i]),na.rm = T)
      delta <- fxsTot - fxTot
      p1 <- fx1[i]/fxsTot
      p2 <- fx2[i]/fxsTot
      p3 <- fx3[i]/fxsTot
      fx1[i] <- fx1[i] - delta*p1
      fx2[i] <- fx2[i] - delta*p2
      fx3[i] <- fx3[i] - delta*p3
    }else{
      fx1[i] <- 0
      fx2[i] <- 0
      fx3[i] <- 0
    }
  }
  return.dx <- data.frame(x, fx, fx1=fx1, fx2=fx2, fx3=fx3)
  return(return.dx)
}

## Function to combine three components
PCLM_eta <- function(x,x1,x2,x3,mx1,mx2,mx3){
  m <- length(x)
  mx1A <- mx2A <- mx3A <- rep(NA,m)
  mx1A[x %in% x1] <- mx1
  mx2A[x %in% x2] <- mx2
  mx3A[x %in% x3] <- mx3
  mx3C <- cbind(mx1A,mx2A,mx3A)
  mxSTAD <- apply(mx3C,1,sum,na.rm=T)
  return(mxSTAD)
}