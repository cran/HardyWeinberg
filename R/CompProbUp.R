CompProbUp <- function(nAA,nBB,EnAB,prob,MaxHet,vec=NULL) {
  prob <-  prob*4*nAA*nBB/((EnAB+2)*(EnAB+1))
  nvec <- c(vec, prob)
  
  while (EnAB < MaxHet-2) {
    nAA <- nAA - 1
    nBB <- nBB - 1
    EnAB <- EnAB + 2
    prob <-  prob*4*nAA*nBB/((EnAB+2)*(EnAB+1))
    nvec <- c(nvec ,prob)
  }
  return(nvec)
}