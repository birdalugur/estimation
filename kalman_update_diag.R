kalman_update_diag <- function(A, C, Q, R, y, x, V, ...){
  # Adapted from programs by Zoubin Ghahramani and Geoffrey E. Hinton, available at http://www.gatsby.ucl.ac.uk/ zoubin, 1996.
  # KALMAN_UPDATE Do a one step update of the Kalman filter
  # [xnew, Vnew, loglik] = kalman_update_diag(A, C, Q, R, y, x, V, ...)
  #
  # INPUTS:
  # A - the system matrix
  # C - the observation matrix 
  # Q - the system covariance 
  # R - the observation covariance
  # y(:)   - the observation at time t
  # x(:) - E[X | y(:, 1:t-1)] prior mean
  # V(:,:) - Cov[X | y(:, 1:t-1)] prior covariance
  #
  # OPTIONAL INPUTS (string/value pairs [default in brackets])
  # 'initial' - 1 means x and V are taken as initial conditions (so A and Q are ignored) [0]
  # 'u'     - u(:) the control signal at time t [ [] ]
  # 'B'     - the input regression matrix
  #
  # OUTPUTS (where X is the hidden state being estimated)
  #  xnew(:) =   E[ X | y(:, 1:t) ] 
  #  Vnew(:,:) = Var[ X(t) | y(:, 1:t) ]
  #  VVnew(:,:) = Cov[ X(t), X(t-1) | y(:, 1:t) ]
  #  loglik = log P(y(:,t) | y(:,1:t-1)) log-likelihood of innovatio#
  # set default params
  
  u <- c()
  B <- c()
  initial <- 0
  
  arguments <- list(...)
  
  for (i in seq(1,length(arguments),by=2)){
    switch (
      arguments[[i]],
      'u' = { u <- arguments[[i+1]] },
      'B' = { B <- arguments[[i+1]] },
      'initial' = { initial <- arguments[[i+1]] },
      { stop('unrecognized argument ',arguments[[i]]) } 
    )
  }
    
  #  xpred(:) = E[X_t+1 | y(:, 1:t)]
  #  Vpred(:,:) = Cov[X_t+1 | y(:, 1:t)]
    
  
  
  if (initial){
    if (isempty(u)){
      xpred <- x
    }
    else { 
      xpred <- x + (B %*% u)
    }
    Vpred <- V
  }
  else{
    if (isempty(u)){
      xpred <- A %*% x
    }
    else{
      xpred <- (A %*% x) + (B %*% u)
    }
    Vpred <- (A %*% V %*% t(A)) + Q
  }
  
  e <- y - C %*% xpred  # error (innovation)
  n <- max(size(e))
  ss <- max(size(A))
  
  d <- size(e,1)
  
  S <- (C %*% Vpred %*% t(C)) + R
  
  GG <- t(C) %*% diag(1 / diag(R)) %*% C
  
  Sinv <- diag(1/diag(R)) - diag(1/diag(R)) %*% C %*% pinv(diag(ss)+(Vpred %*% GG)) %*% Vpred %*% t(C) %*% diag(1/diag(R)) # works only with R diagonal
    
  # Sinv = inv(S);
  
  detS <- prod(diag(R)) %*% det(diag(ss) + Vpred %*% GG)
  
  denom <- (2*pi)^(d/2)*sqrt(abs(detS))
  
  mahal <- apply(t(e) %*% Sinv%*%e,1, sum)
  
  loglik = -0.5%*%mahal - log(denom)
  
  K <- Vpred %*% t(C) %*% Sinv; # Kalman gain matrix

  # If there is no observation vector, set K = zeros(ss).
  xnew <- xpred + K %*% e;              # csi_est(t\t) formula 13.6. 5    
  Vnew <- (diag(ss) - K%*%C)%*%Vpred;    # P(t\t) formula 13.2.16 hamilton
  VVnew <- (diag(ss) - K%*%C)%*%A%*%V;
  
  return(list(xnew=xnew, Vnew=Vnew, loglik=loglik, VVnew=VVnew))
}