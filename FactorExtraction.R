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

source('center.R')
source('kalman_filter_diag.R')
source('kalman_smoother_diag.R')
source('kalman_update_diag.R')
source('ricSW.R')
source('smooth_update.R')




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
  xsmooth <- ksd_result$xsmooth
  Vsmooth <- ksd_result$Vsmooth
  VVsmooth <- ksd_result$VVsmooth
  loglik <- ksd_result$loglik
  # xsmooth = E(F_t)
  # Vsmooth = VAR(F_t)
  
  VF <- Vsmooth
  ind <- size(VF,3)
  F <-  t(xsmooth)
  
  return(list(F=F,VF=VF,A=A,C=C,Q=Q,R=R,initx=initx,initV=initV,ss=ss,MM=MM))
}