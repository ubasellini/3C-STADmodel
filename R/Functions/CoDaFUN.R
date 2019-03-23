

#Smooth from age 80 to 120 with the Kannisto model and remove 0 values at young ages
#Dx = matrix of Death counts by year (row) and age (column)
#Ex = matrix of Exposure to risk by year (row) and age (column)
#mx = matrix of death rates by year (row) and age (column)
fit120<- function(Dx, Ex, mx){
  
  #poisson distribution + kannisto at age 80+
  
  dens<- Dx[,ar.old1]
  Exp<- Ex[,ar.old1]
  age<-c(0:(length(ar.old1)-1))
  
  kan<- function(par,t, D, E){
    a <- par[1]
    b <- par[2]
    haz<- (a*exp(b*t))/(1+(a*exp(b*t)))
    loglike <- (D*log(haz))-(haz*E)
    sumll <- sum(loglike)
    return(sumll)
  }
  
  parameters<-matrix(0, nrow(dens), 2)
  for(i in 1:nrow(parameters)){
    parameters[i,] <- optim(par=c(0.02, 0.1), kan, t=age, D=dens[i,], E=Exp[i,],
                            control = list(fnscale = -1))$par
  }
  
  
  x<-c(0:(length(ar.old2)-1))
  mx.fit<-matrix(0, nrow(dens), length(x))
  for(i in 1:nrow(mx.fit)){
    mx.fit[i,]<-(parameters[i,1]*exp(parameters[i,2]*x))/(1+(parameters[i,1]*exp(parameters[i,2]*x)))
  }
  
  mx80<-mx.fit
  
  
  #0 values
  Dx2<-Dx[,ar.young]
  prop.Dx<- Dx2/rowSums(Dx2)
  prop.Dx[is.na(prop.Dx)]<-0
  for(i in 1:nrow(prop.Dx)){
    for(j in 1:ncol(prop.Dx)){
      if(prop.Dx[i,j]==0){
        prop.Dx[i,j]<- (min(Dx2[Dx2>0], na.rm=T)/2)/sum(Dx2[i,], na.rm=T)
      }
    }
  }
  prop.Dx2<- prop.Dx/rowSums(prop.Dx)
  Dx2<- prop.Dx2*rowSums(Dx[,ar.young])
  mx.young<- Dx2/Ex[,ar.young]
  
  
  mx.adjust<- cbind(mx.young, mx80 )#, mx.fit[,32:41])
  age<-c(ar.young-1, ar.old2-1)
  colnames(mx.adjust)<-age
  rownames(mx.adjust)<-rownames(Ex)
  
  return(mx.adjust)
}


# Life Table starting from dx
#dx = matrix of life table deaths by year (row) and age (column)
LifeT_dx<-function(dx, radix, a, n){
  dx<-dx*radix
  lx<-matrix(radix, nrow(dx), ncol(dx))
  Lx<-matrix(NA, nrow(dx), ncol(dx))
  Tx<-matrix(NA, nrow(dx), ncol(dx))
  for(j in 1:nrow(dx)){
    
    lx[j,2:ncol(dx)]<- radix-cumsum(dx[j,1:(ncol(dx)-1)])
    
    Lx[j,]<- (lx[j,]*n) - (dx[j,]*(n-a))
    
    for(i in 1:ncol(dx)){
      Tx[j,i]<-sum(Lx[j, i:ncol(Lx)], na.rm=T)
    }
    Tx[Tx==0]<-NA
  }
  mx<- dx/Lx
  ex<-Tx/lx
  colnames(lx)<-colnames(Lx)<-colnames(ex)<-colnames(mx)<-colnames(dx)
  rownames(lx)<-rownames(Lx)<-rownames(ex)<-rownames(mx)<-rownames(dx)
  return(list(mx=mx, dx=dx, lx=lx, Lx=Lx, ex=ex))
}


# Life Table starting from mx 
#mx = matrix of death rates by year (row) and age (column)
LifeT_mx<-function(mx, radix, a, n){
  qx<-matrix(NA, nrow(mx), ncol(mx))
  lx<-matrix(NA, nrow(mx), ncol(mx))
  Lx<-matrix(NA, nrow(mx), ncol(mx))
  Tx<-matrix(NA, nrow(mx), ncol(mx))
  for(j in 1:nrow(mx)){
    
    qx[j,]<- mx[j,] * n / ( 1 + ( n - a ) * mx[j,])
    qx[,ncol(qx)]<-1
    qx[qx>1]<-1
    
    px<-1-qx
    
    lx[j,]<- radix* c(1, cumprod(px[j,])[1:(length(lx[j,])-1)])
    
    dx<-lx*qx
    
    Lx[j,]<- (lx[j,]*n) - (dx[j,]*(n-a))
    
    for(i in 1:ncol(mx)){
      Tx[j,i]<-sum(Lx[j,i:ncol(Lx)], na.rm=T)
    }
    Tx[Tx==0]<-NA
  }
  ex<-Tx/lx
  colnames(qx)<- colnames(dx)<-colnames(lx)<-colnames(Lx)<-colnames(Tx)<-colnames(ex)<-colnames(mx)
  rownames(qx)<- rownames(dx)<-rownames(lx)<-rownames(Lx)<-rownames(Tx)<-rownames(ex)<-rownames(mx)
  return(list(mx=mx, dx=dx, lx=lx, Lx=Lx, Tx=Tx, ex=ex))
}

# ThreeDClose function: closing procedure for a 3-dimentional array 
ThreeDClose<-function(boot){
  nrow<-dim(boot)[1]
  ncol<-dim(boot)[2]
  nslice<-dim(boot)[3]
  
  Bclose<-boot
  for(i in 1:nslice){
    Bclose[,,i]<- acomp(boot[,,i])
  }
  return(Bclose)
}

# dx=LTcountry$dx
#Forecast with CoDa
#dx = matrix of life table deaths by year (row) and age (column)
#t = Number of year to be forecast 
CoDa<- function(dx, t, sex){
  
  nrow.dx <- nrow(dx)
  
  #data close
  close.dx<-acomp(dx)
  
  #geometric mean
  ax<- geometricmeanCol(close.dx)
  ax<-ax/sum(ax)
  
  #centering
  dx.cent<-close.dx-acomp(ax)
  
  #clr transformation
  clr.cent<- clr(dx.cent)
  
  # SVD: bx and kt
  par<- svd(clr.cent, nu=1, nv=1)
  U <- par$u
  V <- par$v
  S <- diag(par$d)
  
  bx<- V[,1]
  kt<- S[1,1]*U[,1]
  
  
  # forecast kt
  TS <- Arima(kt, order=order.coda ,include.drift=T, method="ML")
  fcst <- forecast(TS, h=t)
  kt.fit <- fcst$fitted
  kt.for <- fcst$mean
  
  variability<- cumsum((par$d)^2/sum((par$d)^2))
  
  #correcting kt trend due to the break created by the MA component. Do not used this approach for other ARIMA order than (0,1,1).
  pert<- (kt.for[2]-kt.for[1]) - (kt.for[1] - kt[length(kt)] )
  kt.for.adj<-kt.for+pert
  
  #projections
  tbx <- t(bx)
  clr.proj.fit<- matrix(kt.fit,(nrow.dx),1) %*% tbx
  clr.proj.for<- matrix(c(kt, kt.for.adj),nrow.dx+t,1) %*% tbx
  
  
  #bootstrap
  bootstrap <- array(NA, c(t, ncol(clr.proj.for), n.kt*n.error))
  error<- clr.proj.fit-clr.cent
  
  kseq <- jseq <- 1:n.kt # kseq defines slices in 3D; jseq defines columns in 2D
  
  for(i in 1:n.error){
    par2<- svd(sample(error, replace=T) + clr.proj.fit, nu=1, nv=1)
    tv1 <- t(par2$v[,1])
    
    test <- Arima(diag(par2$d)[1,1]*par2$u[,1], order=order.coda, include.drift=T, method="ML")
    # fill a matrix with forecast simulations as columns
    boot.kt <- replicate(n.kt, simulate(test, nsim=t, bootstrap=T), simplify = "array")
    
    # simplify = "array" gives 3D result; sapply defaults to 2D
    bootstrap[ , , kseq] <-
      sapply(jseq, function(j) (boot.kt[ ,j] %*% tv1), simplify = "array")
    
    kseq <- kseq + n.kt # define next set of slices
    
  }
  
  
  #Inverse clr transformation
  BK.proj.fit<- clrInv(clr.proj.fit)
  BK.proj.for<- clrInv(clr.proj.for)
  BK.proj.boot <- ThreeDClose(exp(bootstrap))
  rm(bootstrap)
  
  
  #Add back geometric mean
  proj.fit<- BK.proj.fit+acomp(ax)
  proj.for<- BK.proj.for+acomp(ax)
  
  proj.boot<-sweep(BK.proj.boot, c(2,1), ax, FUN="*")
  close.proj.boot <- ThreeDClose(proj.boot)
  rm(proj.boot)
  
  #jump-off correction
  rows <- nrow.dx:nrow(proj.for)
  pert<- acomp(dx[nrow.dx,])-acomp(proj.for[nrow.dx,])
  forecast.dx<- acomp(proj.for[rows,]) +pert
  forefit.dx<- proj.for
  forefit.dx[rows,]<-forecast.dx
  
  
  #jump-off for bootstrap
  median<-apply(close.proj.boot, c(1,2), function(close.proj.boot) quantile(close.proj.boot, prob=0.5, type=8))
  median.close<- median/rowSums(median)
  
  pert.boot<- forefit.dx[(nrow.dx+1),]/median.close[1,]
  pert.boot.close<-pert.boot/sum(pert.boot)
  forecast.boot<- sweep(close.proj.boot, c(2,1), pert.boot.close, FUN="*")
  rm(close.proj.boot)
  close.forecast.boot <- ThreeDClose(forecast.boot)
  rm(forecast.boot)
  
  
  #Life table calculation
  LT<-LifeT_dx(unclass(forefit.dx), radix, a, n)
  lx<-LT$lx
  mx<-LT$mx
  ex<-LT$ex
  
  ## LE 
  ex.int <- apply(close.forecast.boot,3,
                  function(close.forecast.boot) LifeT_dx(close.forecast.boot, radix, a, n)$ex)
  
  dim(ex.int) <- dim(close.forecast.boot)
  
  
  int95_ex <- apply(ex.int, c(1,2),
                    function(ex.int) quantile(ex.int, prob=prob, type=8))
  
  ## MX 
  mx.int <- apply(close.forecast.boot,3,
                  function(close.forecast.boot) LifeT_dx(close.forecast.boot, radix, a, n)$mx)
  
  dim(mx.int) <- dim(close.forecast.boot)
  
  
  int95_mx <- apply(mx.int, c(1,2),
                    function(mx.int) quantile(mx.int, prob=prob, type=8))
  
  ## GINI
  gx.int <- matrix(NA,dim(mx.int)[1],dim(mx.int)[3])
  for (jj in 1:dim(mx.int)[3]){
    gx.int[,jj] <- apply(mx.int[,,jj],1,GINI_func,ages=age2,sex=substr(sex,1,1))
    
  }
  
  int95_gx <- apply(gx.int, 1,
                    function(gx.int) quantile(gx.int, prob=prob, type=8))
  
  
  #Return list
  return(list(mx=mx, dx=forefit.dx, ex=ex,  
              ex.lw95=int95_ex[1,,1], ex.lw80=int95_ex[2,,1],  ex.up95=int95_ex[5,,1], ex.up80=int95_ex[4,,1],
              ex.med=int95_ex[3,,1], 
              gx.lw95=int95_gx[1,], gx.lw80=int95_gx[2,],  gx.up95=int95_gx[5,], gx.up80=int95_gx[4,],
              gx.med=int95_gx[3,],
              mx.lw95=int95_mx[1,,], mx.lw80=int95_mx[2,,],  mx.up95=int95_mx[5,,], mx.up80=int95_mx[4,,],
              mx.med=int95_mx[3,,],
              VarExp=variability, ax=ax, bx=bx, kt=kt,
              BK=BK.proj.for,BKsim=BK.proj.boot)) #BK and BKsim for the average have to be extracted and included in the CoDaCh function
}
