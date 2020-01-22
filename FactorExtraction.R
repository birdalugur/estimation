FactorExtraction <- function(x,q,r,p,A,C,Q,R,initX,initV,ss,MM){
  "
  x <- matrix(1:12, nrow = 3, ncol = 4)
  x : Matrix "
  
  T <- dim(x)[1]
  N <- dim(x)[2]
  
  das <- colSums(is.na(x))
  
  m <- max(das)
  
  if (nargs()<5){
    #statement1
  }
  else{
    #statement2
  }
  
}