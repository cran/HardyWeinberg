`HWLratio` <-
function(X,verbose=FALSE)
{
  nAA <- X[1]
  nAB <- X[2]
  nBB <- X[3]
  n <- sum(X)
  nA <- 2*nAA+nAB
  nB <- 2*nBB+nAB
  lnLambda <- n*log(n) + nAB*log(2) + nA*log(nA) + nB*log(nB) - 2*n*log(2*n) - nAA*log(nAA) - nAB*log(nAB) - nBB*log(nBB)
  Lambda <- exp(lnLambda)
  G2 <- -2*lnLambda
  pval <- pchisq(G2,1,lower.tail=FALSE)
  if(verbose) cat("G2 =",G2,"p-value =",pval,"\n")
  return(list(Lambda=Lambda,G2=G2,pval=pval))
}

