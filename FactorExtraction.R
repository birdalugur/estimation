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
}


center <- function(x){
  T <- dim(x)[1]
  N <- dim(x)[2]
  xc <- x-(matrix(1, nrow = T, ncol = N) %*% diag(apply(z, 2, sum)/T))
  return(xc)
}