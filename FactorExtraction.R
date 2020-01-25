library(pracma)
library(Matrix)

FactorExtraction <- function(x,q,r,p,A,C,Q,R,initX,initV,ss,MM){
  "
  x <- matrix(1:12, nrow = 3, ncol = 4)
  x : Matrix
  T : number of columns
  N : number of rows
  das : number of nan's in columns
  m : maximum value in das"
  
  T <- dim(x)[1]
  N <- dim(x)[2]
  
  das <- colSums(is.na(x))
  
  m <- max(das)
  
  if (nargs()<5){
    selected_num_row <- T-m
    z <- x[1:selected_num_row,]
    ss <- apply(z, 2, sd)
    MM <- apply(z, 2, mean)
    s <- matrix(1, nrow = T, ncol = length(ss)) %*% diag(ss)
    M <- matrix(1, nrow = T, ncol = length(ss)) %*% diag(MM)
    x = (x - M)/s
    z <- x[1:selected_num_row,]
    
    #[A, C, Q, R, initx, initV] = ricSW(z,q,r,p);
  }
  else{
    s <- matrix(1, nrow = T, ncol = length(ss)) %*% diag(ss)
    M <- matrix(1, nrow = T, ncol = length(ss)) %*% diag(MM)
    x = (x - M)/s
    z <- x[1:selected_num_row,]
  }
  
}


ricSW <- function(x,q,r,p){
  Mx=apply(x, 2, mean)
  Wx=diag(apply(x, 2, sd))
  x=center(x)%*%inv(Wx)
  
  OPTS.disp = 0
    
  T <- dim(x)[1]
  N <- dim(x)[2]
  
  if (r < q) {
    message('q has to be less or equal to r')
  }
  
  nlag <- p-1
  A_temp <- matrix(0L, nrow = r, ncol = r*p) 
  I <- diag(r*p)
  end_I <- dim(I)[1]-r
  if (end_I != 0){
    A <- rbind(t(A_temp),I[1:end_I-r,])
  }
  else{
    A <- rbind(t(A_temp),I[0,])
  }
  
  Q <- matrix(0L, nrow = r*p, ncol = r*p) 
  Q[1:r,1:r] <- diag(r)
  OPTS.disp = 0
  
  "NOT -------------d'nin köşegen matris olması gerekiyor"
  result_eigs <- eigs(cov(A),k=2,which = "LM")	# computes eigenvalues and eigenvectors of the var-covariance 
  d <- result_eigs$values
  v <- result_eigs$vectors
  
  # matrix of the data, x.
  # d is a rxr diagonal matrix with the 10 largest eigenvalues on the diagonal. 
  # v is a nxr matrix of the eigenvectors that corresponds to the eigenvalues.
}


center <- function(x){
  T <- dim(x)[1]
  N <- dim(x)[2]
  xc <- x-(matrix(1, nrow = T, ncol = N) %*% diag(apply(x, 2, sum)/T))
  return(xc)
}



