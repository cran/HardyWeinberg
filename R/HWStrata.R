HWStrata <- function(X,verbose=TRUE) {
  X  <- X[,names(order.auto(X[1,]))]
  K  <- nrow(X) # number of strata
  
  nAA <- X[,1]
  nAB <- X[,2]
  nBB <- X[,3]
  
  nk <- rowSums(X)
  pk <- af(X)
  qk <- 1 - pk
  
  num.1   <- nAB*(nAB - 1) - 4*nAA*nBB
  den.1   <- 2*(2*nk-1)
  ratio.1 <- num.1/den.1
  NUM <- (sum(ratio.1))^2
 
  p2q2  <- ((8*nk*nk*nk)/((2*nk-1)*(2*nk-2)*(2*nk-3)))*pk*pk*qk*qk
  p2q2  <- p2q2 - ((2*nk)/((2*nk-2)*(2*nk-3)))*pk*qk
  num.2 <- 4*nk*(nk-1)*p2q2
  den.2 <- den.1
  ratio.2 <- num.2/den.2
  DEN <- sum(ratio.2)
  
  T2 <- NUM/DEN
  
  pval <- pchisq(T2,df=1,lower.tail = FALSE)
  
  if(verbose) {
    cat("Olson's asymptotic test for HWE across strata for a single biallelic marker.\n")
    cat("T2 = ",T2,"p-value = ",pval,"\n")
  }
  
  out <- list(T2=T2,pval=pval)
}
