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
      arguments[[i]],
      'model' = { model <- arguments[[i+1]]  }, 
      'u' = { u <- arguments[[i+1]] },
      'B' = { B <- arguments[[i+1]] },
      { stop('unrecognized argument ',arguments[[i]]) } 
    )
    
  }
  
  xsmooth <- matrix(0L, nrow = ss, ncol = T)
  Vsmooth <- (array(c(0),dim = c(ss,ss,T)))
  VVsmooth <- (array(c(0),dim = c(ss,ss,T)))
  
  # Forward pass
  kfd_result <- kalman_filter_diag(y, A, C, Q, R, init_x, init_V, 'model', model, 'u', u, 'B', B)  #burda kaldık
  xfilt  <- kfd_result$x
  Vfilt  <- kfd_result$V
  VVfilt <- kfd_result$VV
  loglik <- kfd_result$loglik
  
  'KONTROL EDİLECEK !!!!!!!!!!!!!!!!!!!!!!!!!!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<< LINE 171'
  # Backward pass
  xsmooth[,T] = xfilt[,T]
  Vsmooth[,,T] = Vfilt[,,T]
  #VVsmooth[,,T] = Vfilt[,,T]
  #browser()
  for (t in T:1) {
    m <- model[t+1]
    print(A[,,m])
    if (isempty(B)){
      result_s_update = smooth_update(
        xsmooth[,t+1,drop=FALSE],
        Vsmooth[,,t+1],
        xfilt[,t],
        Vfilt[,,t],
        Vfilt[,,t+1],
        VVfilt[,,t+1],
        A[,,m],
        Q[,,m],
        c(),
        c()
      )
      xsmooth[,t,drop=FALSE] <- result_s_update$xsmooth
      Vsmooth[,,t] <- result_s_update$Vsmooth
      VVsmooth[,,t] <- result_s_update$VVsmooth_future
    }
    
    else{
      result_s_update <- smooth_update(
        xsmooth[,t+1,drop=FALSE],
        Vsmooth[,,t+1],
        xfilt[,t],
        Vfilt[,,t],
        Vfilt[,,t+1],
        VVfilt[,,t+1],
        A[,,m],
        Q[,,m],
        B[,,m],
        u[,,t+1]
      )
      xsmooth[,t,drop=FALSE] <- result_s_update$xsmooth
      Vsmooth[,,t] <- result_s_update$Vsmooth
      VVsmooth[,,t] <- result_s_update$VVsmooth_future
      
    }
  }
  
  VVsmooth[,,] = matrix(0L,ss,ss)
  
  return(list(xsmooth=xsmooth, Vsmooth=Vsmooth, VVsmooth=VVsmooth, loglik=loglik))
  
}