library(pracma)
library(Matrix)
library(RSpectra)

# NOT ----------
# eigs fonksiyonları kontrol edilmeli

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
  d <- diag(length(result_eigs$values))*result_eigs$values
  v <- result_eigs$vectors
  "v nin 1. sütununun işareti yanlış -> v[,1]=v[,1]*-1 "
  v[,1]=v[,1]*-1
  
  # matrix of the data, x.
  # d is a rxr diagonal matrix with the 10 largest eigenvalues on the diagonal. 
  # v is a nxr matrix of the eigenvectors that corresponds to the eigenvalues.
  
  F <- x%*%v
  
  ' cov fonksiyonu aynı girdiye farklı çıktılar verebilir(matlab cov dan farklı) !'
  R <- diag(diag(cov((x-x%*%v%*%t(v))))) #Estimate of the covariance matrix of the idiosincratic component
  # REMARK: x*v*v' is the projection of x over the principal components (F=x*v)

  if (p>0) { 
    #ESTIMATE the AUTOregressive model for the Factors: run the var F(t) = A_1*F(t-1)+...+A_p*F(t-1) + e(t);
    z = F
    "Matris satırlarının sayısı eşleşmeli ! ! ! !"
    ############################################# HATALI ##
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
    for (kk in 1:nlag){  'for döngüsü kontrol edilmeli'
      Z <- cbind(Z,z[(nlag-kk+1):(size(z)[1]-kk),]) # stacked regressors (lagged SPC)
      
    }
    initx <- t(t(Z[1,]))
    initV <- matrix((pinv(diag(size(kron(A,A),1))-kron(A,A)) %*% (matrix(as.vector(Q),  ncol = 1))),r*p,r*p) # initV = cov(Z); %eye(r*(nlag+1))
  }
  else{
    initx <- c()
    initv <- c()
  }
  C <- cbind(v,matrix(0L, nrow = N, ncol = r*nlag))
}


center <- function(x){
  T <- dim(x)[1]
  N <- dim(x)[2]
  xc <- x-(matrix(1, nrow = T, ncol = N) %*% diag(apply(x, 2, sum)/T))
  return(xc)
}



