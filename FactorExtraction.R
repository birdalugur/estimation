setwd("/home/ugur/r_projects/estimation/")
library(pracma)
library(Matrix)
library(RSpectra)
library(rmatio)


data<-read.mat('data.mat')
x <- data$X
p <- data$p
q <- data$q
r <- data$r

FactorExtraction <- function(x,q,r,p,A,C,Q,R,initX,initV,ss,MM){
  "
  x <- matrix(1:12, nrow = 3, ncol = 4)
  x : Matrix
  T : number of columns
  N : number of rows
  das : number of nan's in columns
  m : maximum value in das"
  
  # function [F,VF,A,C,Q,R,initx,initV,ss,MM] = FactorExtraction(x,q,r,p,A,C,Q,R,initx,initV,ss,MM);
  # extract common factors from vector of time series possibly unbalanced
  # at the end of the sample, (NaN for missing observations)
  
  #The model
  #x_t = C F_t + \xi_t
  #F_t = AF_{t-1} + B u_t
  #R = E(\xi_t \xi_t')
  # Q = BB'
  # u_t ~ WN(0,I_q)
  # initx = F_0
  # initV = E(F_0 F_0')
  # ss: std(x) 
  # MM: mean(x)
  
  # q: dynamic rank
  # r: static rank (r>=q)
  # p: ar order of the state vector (default p=1)
  
  # F : estimated factors
  # VF: estimation variance for the common factors
  
  T <- dim(x)[1]  # dimension of the panel
  N <- dim(x)[2]
  
  # Construct the balanced panel z from the original panel x
  # NOTES: sum(isnan(x)) computes the number of NaNs in each column 
  # of x and stores that number in a cell in a row vector, das.
  
  das <- colSums(is.na(x))
  
  m <- max(das)
  
  if (nargs()<5){
    # Estimate the parameters, if they are not provided, by simple regrssion on
    # Principal components estimates of the common factors (based on the balanced part of the panel)
    
    
    selected_num_row <- T-m
    z <- x[1:selected_num_row,] # so z is the matrix with # of rows = T-m (all rows with any NaNs are excluded)
    ss <- apply(z, 2, sd) # computes stdev of each column of data.
    MM <- apply(z, 2, mean)
    
    # STEP:  Standardize the panel
    s <- matrix(1, nrow = T, ncol = length(ss)) %*% diag(ss)
    M <- matrix(1, nrow = T, ncol = length(ss)) %*% diag(MM)
    x = (x - M)/s
    z <- x[1:selected_num_row,]
    
    
    result_ricsw <- ricSW(z,q,r,p)
    A <- result_ricsw$A
    C <- result_ricsw$C
    Q <- result_ricsw$Q
    R <- result_ricsw$R
    initx <- result_ricsw$initx
    initV <- result_ricsw$initV
  }
  else{ #if the parameters are given, just standardize the variables with the mean and std computed over the balanced panel.
    s <- matrix(1, nrow = T, ncol = length(ss)) %*% diag(ss)
    M <- matrix(1, nrow = T, ncol = length(ss)) %*% diag(MM)
    x = (x - M)/s
    z <- x[1:selected_num_row,]
  }
  
  # The signal extraction in presence of missing data is performed by
  # using a time varying Kalman filter in which missing data are assigned an
  # extremely large variance of the noise in idiosyncratic component.
  
  # Define the parameters of the time varying state space model... time is
  # on the 3rd dimension
  
  AA <- (array(c(A),dim = c(dim(A),T)))
  QQ <- (array(c(Q),dim = c(dim(Q),T)))
  CC <- (array(c(C),dim = c(dim(C),T)))
  RR <- (array(c(R),dim = c(dim(R),T)))
  
  
  for (jt in 1:T){
    miss <- is.nan(x[jt,])
    Rtemp <- matrix(diag(R),ncol = 1)
    Rtemp[miss] = 1e+32
    RR[,,jt] =  diag(c(Rtemp))
  }
  
  
  xx<- x
  xx[is.nan(x)] = 0 # missing data are assigned an arbitrary value...
  
  #  Run the kalman smoother on the time varying state space model
  ksd_result <- kalman_smoother_diag(t(xx),AA, CC, QQ, RR, initx, initV,'model',1:T)
  xsmooth <- 1
  Vsmooth <- 1
  VVsmooth <- 1
  loglik <- 1
  # xsmooth = E(F_t)
  # Vsmooth = VAR(F_t)
  
  VF <- Vsmooth
  ind <- size(VF,3)
  F <-  t(xsmooth)
}

kalman_smoother_diag <- function(y, A, C, Q, R, init_x, init_V, ...){
  # Adapted from programs by Zoubin Ghahramani and Geoffrey E. Hinton, available at http://www.gatsby.ucl.ac.uk/ zoubin, 1996.
  # Kalman/RTS smoother.
  # [xsmooth, Vsmooth, VVsmooth, loglik] = kalman_smoother_diag(y, A, C, Q, R, init_x, init_V, ...)
  #
  # The inputs are the same as for kalman_filter.
  # The outputs are almost the same, except we condition on y(:, 1:T) (and u(:, 1:T) if specified),
  # instead of on y(:, 1:t).
  
  os <- size(y)[1]
  T <- size(y)[2]
  ss <- size(A,1)
  
  # set default params
  
  model <- matrix(1, nrow = 1, ncol = T)
  
  u <- c()
  B <- c()
  
  arguments <- list(...)
  n_arguments <- length(arguments)
  
  for (i in seq(1,n_arguments,by=2)){
    switch (
      arguments[1],
            'model' = { model <- arguments[i+1] },
            'u' = { u <- arguments[i+1] },
            'B' = { B <- arguments[i+1] },
            { stop('unrecognized argument ',arguments[i]) } 
      )
    
  }
  
  xsmooth <- matrix(0L, nrow = ss, ncol = T)
  Vsmooth <- (array(c(0),dim = c(2,2,T)))
  VVsmooth <- (array(c(0),dim = c(2,2,T)))
  
  # Forward pass
  kfd_result <- kalman_filter_diag(y, A, C, Q, R, init_x, init_V, 'model', model, 'u', u, 'B', B)
  xfilt  <- kfd_result$xfilt
  Vfilt  <- kfd_result$Vfilt
  VVfilt <- kfd_result$VVfilt
  loglik <- kfd_result$loglik
  
  'KONTROL EDİLECEK !!!!!!!!!!!!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<< LINE 171'
  # Backward pass
  xsmooth[,T] = xfilt[,T]
  Vsmooth[,,T] = Vfilt[,,T]
  #VVsmooth[,,T] = Vfilt[,,T]
  
  for (t in T:1) {
    
  }
  
}

ricSW <- function(x,q,r,p){
  # %ricSW(z,q,r,p);
  # Computes the parameters of the factor models 
  # REMARK: the parameters C and R refer to the standardized variables.
  
  Mx=apply(x, 2, mean) # Mean
  Wx=diag(apply(x, 2, sd)) # Standard deviation
  x=center(x)%*%inv(Wx) # Standardize
  
  OPTS.disp = 0
    
  T <- dim(x)[1]  # size of the database
  N <- dim(x)[2]
  
  if (r < q) {
    message('q has to be less or equal to r')
  }
  
  nlag <- p-1 # p=1, so nlag = 0.
  
  # Define some preliminary quantity that are necessary to writhe the VAR in companion form
  A_temp <- matrix(0L, nrow = r, ncol = r*p)  # a zero matrix,
  I <- diag(r*p) # identity matrix,
  
  # NOTE: if p=1, then I(1:end-r,1:end) is empty. In this case, MATLAB reads A as equal to A_temp.
  end_I <- dim(I)[1]-r
  if (end_I != 0){
    A <- rbind(t(A_temp),I[1:end_I-r,])
  }
  else{
    A <- rbind(t(A_temp),I[0,])
  }
  
  Q <- matrix(0L, nrow = r*p, ncol = r*p)  #a zero matrix, 10x10.
  Q[1:r,1:r] <- diag(r) #identity of size=10. 
  OPTS.disp = 0
  
  
  result_eigs <- eigs(cov(x),k=r,which = "LM")	# computes eigenvalues and eigenvectors of the var-covariance 
  # matrix of the data, x.
  # d is a rxr diagonal matrix with the 10 largest eigenvalues on the diagonal. 
  # v is a nxr matrix of the eigenvectors that corresponds to the eigenvalues.
  d <- diag(length(result_eigs$values))*result_eigs$values
  v <- result_eigs$vectors
  v[,1]=v[,1]*-1
  
  
  
  F <- x%*%v  # PC estimates of the common factors
  
  ' cov fonksiyonu aynı girdiye farklı çıktılar verebilir(matlab cov dan farklı) !'
  R <- diag(diag(cov((x-x%*%v%*%t(v))))) #Estimate of the covariance matrix of the idiosincratic component
  # REMARK: x*v*v' is the projection of x over the principal components (F=x*v)

  if (p>0) { 
    #ESTIMATE the AUTOregressive model for the Factors: run the var F(t) = A_1*F(t-1)+...+A_p*F(t-1) + e(t);
    z = F
    #Z<-matrix(, nrow = size(z)[1], ncol = 0)
    Z <- c()
    for (kk in 1:1){
      Z <- cbind(Z,z[(p-kk+1):(size(z)[1]-kk),]) # stacked regressors (lagged SPC)
    }
    ##############################################
    z<-z[(p+1):size(z)[1],]
    A_temp <- (inv(t(Z)%*%Z)%*%t(Z))%*%z #OLS estimator of the VAR transition matrix
    A[1:r,1:r*p] = t(A_temp) 
    
    # Compute Q
    e <- z-Z%*%A_temp # VAR residuals
    H <- cov(e) # VAR covariance matrix
    
    if (r==q){
      # The covariance matrix of the VAR residuals is of full rank
      Q[1:r,1:r] = H
    }
    else{ 'Bu blok kontrol edilmeli '
      # The covariance matrix of the VAR residuals has reduced rank
      res_ed <- eigs(H,k=2,which = "LM") # eigenvalue decomposition
      P <- res_ed$vectors
      M <- res_ed$values
      M <- diag(length(M))*M
      
      # P<- matrix(c(-0.9530,0.3029,-0.3029,-0.9530), nrow = 2,ncol = 2,byrow = TRUE) 
      # M<-matrix(c(1.7018,0,0,1.1271),nrow = 2,ncol = 2)
      
      P <- e%*%P%*%diag(sign(P[1,]))
      M ^ (-0.5)
      
      "matrix power- (+) tamsayı olmayan k değerinden dolayı hatalı sonuç vermekte"
      'https://www.rdocumentation.org/packages/expm/versions/0.999-4/topics/matpow'
      # u_orth = e*P*(M^-.5); # extracting the common shocks
      
      
      e_pc <- e %*% P %*% t(P)
      Q[1:r,1:r] = P %*% M %*% t(P)
    }
  }
  
  # Computes the initial conditions for the filter.
  # The common factors are initialized by the PC estimates.
  # Initial variance is set equal to the unconditional variance ofthe common factors. 
  
  if (p > 0){
    z <- F
    Z <- c()
    for (kk in 0:nlag){  
      Z <- cbind(Z,z[(nlag-kk+1):(size(z)[1]-kk),]) # stacked regressors (lagged SPC)
      
    }
    initx <- t(t(Z[1,]))
    initV <- matrix((pinv(diag(size(kron(A,A),1))-kron(A,A)) %*% (matrix(as.vector(Q),  ncol = 1))),r*p,r*p) # initV = cov(Z); %eye(r*(nlag+1))
  }
  else{
    initx <- c()
    initv <- c()
  }
  C <- cbind(v,matrix(0L, nrow = N, ncol = r*nlag)) # Cov(data,factors); recall nlag = 0.
  
  return(list(A=A, C=C, Q=Q, R=R, initx=initx, initV=initV, Mx=Mx, Wx=Wx))
}


center <- function(x){
  T <- dim(x)[1]
  N <- dim(x)[2]
  xc <- x-(matrix(1, nrow = T, ncol = N) %*% diag(apply(x, 2, sum)/T))
  return(xc)
}
