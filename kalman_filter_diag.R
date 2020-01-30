kalman_filter_diag <- function(y, A, C, Q, R, init_x, init_V, ...){
  # Adapted from programs by Zoubin Ghahramani and Geoffrey E. Hinton, available at http://www.gatsby.ucl.ac.uk/ zoubin, 1996.
  # Kalman filter.
  # [x, V, VV, loglik] = kalman_filter_diag(y, A, C, Q, R, init_x, init_V, ...)
  #
  # INPUTS:
  # y(:,t)   - the observation at time t
  # A - the system matrix
  # C - the observation matrix 
  # Q - the system covariance 
  # R - the observation covariance
  # init_x - the initial state (column) vector 
  # init_V - the initial state covariance 
  #
  # OPTIONAL INPUTS (string/value pairs [default in brackets])
  # 'model' - model(t)=m means use params from model m at time t [ones(1,T) ]
  #     In this case, all the above matrices take an additional final dimension,
  #     i.e., A(:,:,m), C(:,:,m), Q(:,:,m), R(:,:,m).
  #     However, init_x and init_V are independent of model(1).
  # 'u'     - u(:,t) the control signal at time t [ [] ]
  # 'B'     - B(:,:,m) the input regression matrix for model m
  #
  # OUTPUTS (where X is the hidden state being estimated)
  # x(:,t) = E[X(:,t) | y(:,1:t)]
  # V(:,:,t) = Cov[X(:,t) | y(:,1:t)]
  # VV(:,:,t) = Cov[X(:,t), X(:,t-1) | y(:,1:t)] t >= 2
  # loglik = sum{t=1}^T log P(y(:,t))
  #
  # If an input signal is specified, we also condition on it:
  # e.g., x(:,t) = E[X(:,t) | y(:,1:t), u(:, 1:t)]
  # If a model sequence is specified, we also condition on it:
  # e.g., x(:,t) = E[X(:,t) | y(:,1:t), u(:, 1:t), m(1:t)]
  
  
  os <- size(y)[1]
  T <- size(y)[2]
  ss <- size(A,1) #size of state space
  
  # set default params
  
  model <- matrix(1, nrow = 1, ncol = T)
  
  u <- c()
  B <- c()
  ndx <- c()
  
  arguments <- list(...)
  n_arguments <- length(arguments)
  
  for (i in seq(1,n_arguments,by=2)){
    switch (
      arguments[[i]],
      'model' = { model <- arguments[[i+1]] },
      'u' = { u <- arguments[[i+1]] },
      'B' = { B <- arguments[[i+1]] },
      'ndx' = { ndx <- arguments[[i+1]] },
      { stop('unrecognized argument ',arguments[[i]]) } 
    )
    
  }
  
  
  x <- matrix(0L, nrow = ss, ncol = T)
  V <- (array(c(0),dim = c(ss,ss,T)))
  VV <- (array(c(0),dim = c(ss,ss,T)))
  
  loglik <- 0
  
  'BURADA m DEĞERİNİ 0 ALGILIYOR ! '
  for (t in 1:T){
    m <- model[t]
    print(m)
    print('-----------------')
    print(A[, ,m])
    if (t == 1){
      # prevx <- init_x[,m]
      # prevV <- init_V[, ,m]
      prevx <- init_x
      prevV <- init_V
      initial <- 1
    }
    else{
      prevx <- x[,t-1,drop=FALSE]     
      prevV <- V[, ,t-1]
      initial <- 0
    }
    
    if (isempty(u)){
      result_kud <- kalman_update_diag(
        A[, ,m],
        C[, ,m],
        Q[, ,m],
        R[, ,m],
        y[, t],
        prevx,
        prevV,
        'initial', 
        initial
      )
      x[,t] <- result_kud$xnew
      V[,,t] <- result_kud$Vnew
      LL <- result_kud$loglik
      VV[,,t] <-  result_kud$VVnew
    }
    else{
      if (isempty(ndx)){
        result_kud <- kalman_update_diag(
          A[, ,m],
          C[, ,m],
          Q[, ,m],
          R[, ,m],
          y[, t],
          prevx,
          prevV,
          'initial', 
          initial,
          'u',
          u[,t],
          'B',
          B[,,m]
        )
        x[,t,drop=FALSE] <- result_kud$xnew
        V[,,t] <- result_kud$Vnew
        LL <- result_kud$loglik
        VV[,,t] <- result_kud$VVnew
      }
      else{
        "ndx köşeli parantez kontrol edilecek <- <<-< -< -<- -< -<- -<- <- -<- -<-<--<-<-<-<-<-<--<-<-<-<--<-<-<-<--<"
        i <- ndx[t]
        # copy over all elements; only some will get updated
        x[,t] <- prevx
        prevP <- inv(prevV)
        prevPsmall <- prevP[i,i]
        prevVsmall <- inv(prevPsmall)
        
        result_kud <- kalman_update_diag(
          A[i,i,m],
          C[,i,m],
          Q[i,i,m],
          R[, ,m],
          y[, t],
          prevx[i],
          prevVsmall,
          'initial', 
          initial,
          'u',
          u[,t],
          'B',
          B[i,,m]
        )
        
        x[i,t] <-result_kud$xnew
        smallV <-result_kud$Vnew
        LL <-result_kud$loglik
        VV[i,i,t] <-result_kud$VVnew
          
          
          
        smallP <- inv(smallV)
        prevP[i,i] <- smallP
        V[,,t] <- inv(prevP)
      }
    }
    loglik <- loglik + LL
  }
  
  
  return(list(x=x, V=V, VV=VV, loglik=loglik))
}