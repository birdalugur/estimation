smooth_update <- function(xsmooth_future, Vsmooth_future,xfilt, Vfilt,  Vfilt_future, VVfilt_future, A, Q, B, u){
  # Adapted from programs by Zoubin Ghahramani and Geoffrey E. Hinton, available at http://www.gatsby.ucl.ac.uk/ zoubin, 1996.
  # One step of the backwards RTS smoothing equations.
  # function [xsmooth, Vsmooth, VVsmooth_future] = smooth_update(xsmooth_future, Vsmooth_future, ...
  #    xfilt, Vfilt,  Vfilt_future, VVfilt_future, A, B, u)
  #
  # INPUTS:
  # xsmooth_future = E[X_t+1|T]
  # Vsmooth_future = Cov[X_t+1|T]
  # xfilt = E[X_t|t]
  # Vfilt = Cov[X_t|t]
  # Vfilt_future = Cov[X_t+1|t+1]
  # VVfilt_future = Cov[X_t+1,X_t|t+1]
  # A = system matrix for time t+1
  # Q = system covariance for time t+1
  # B = input matrix for time t+1 (or [] if none)
  # u = input vector for time t+1 (or [] if none)
  #
  # OUTPUTS:
  # xsmooth = E[X_t|T]
  # Vsmooth = Cov[X_t|T]
  # VVsmooth_future = Cov[X_t+1,X_t|T]#
  # xpred = E[X(t+1) | t]
  
  if (isempty(B)){
    xpred <- A %*% xfilt
  }
  else{
    xpred <- (A %*% xfilt) + B %*% u
  }
  
  Vpred <- A %*% Vfilt %*% t(A) + Q  # Vpred = Cov[X(t+1) | t]
  
  J <- Vfilt %*% t(A) %*% pinv(Vpred) # smoother gain matrix
  
  xsmooth <- xfilt + J %*% (xsmooth_future - xpred)
  
  Vsmooth <- Vfilt + J %*% (Vsmooth_future - Vpred) %*% t(J)
  
  VVsmooth_future <- VVfilt_future + (Vsmooth_future - Vfilt_future) %*% pinv(Vfilt_future) %*% VVfilt_future

  return(list(xsmooth=xsmooth, Vsmooth=Vsmooth, VVsmooth_future=VVsmooth_future))
  
}