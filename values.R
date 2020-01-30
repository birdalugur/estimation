y=t(xx)
A=AA
C=CC
Q=QQ 
R=RR
init_x=initx
init_V=initV
varargin = list('model',1:T)



varargin = list('model', model, 'u', u, 'B', B)
arguments <- varargin


# kalman update diag
A = A[,,m]
C = C[,,m]
Q = Q[,,m]
R = R[,,m]
y = y[,t]
x = prevx
V = prevV
varargin = list('initial', initial,'u',u[,t],'B',B[,,m]) 
