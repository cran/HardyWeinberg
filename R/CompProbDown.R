CompProbDown <- function(nAA,nBB,EnAB,prob,vec=NULL) 
{
  prob <-  prob*EnAB*(EnAB-1)/(4*(nAA+1)*(nBB+1))
  nvec <- c(vec, prob)
  
  while(EnAB > 3) {
    nAA <- nAA + 1
    nBB <- nBB + 1
    EnAB <- EnAB - 2
    prob <-  prob*EnAB*(EnAB-1)/(4*(nAA+1)*(nBB+1))
    nvec <- c(nvec ,prob)
  }
  
  return(nvec)
}