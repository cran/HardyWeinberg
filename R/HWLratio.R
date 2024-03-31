HWLratio <- function (X, verbose = TRUE, x.linked = FALSE) {
if(!x.linked) { # autosomal marker
  if (length(X) != 3 | any(X < 0)) 
    stop("HWLratio: X is not a 3 by 1 non-negative count vector")
  if (any(!is.wholenumber(X))) {
    warning("Genotype counts are not integers, counts will be rounded.")
    X <- round(X, digits = 0)
  }
  X <- order.auto(X)
  n <- sum(X)
  nA <- 2 * X[1] + X[2]
  nB <- 2 * X[3] + X[2]
  nxxvec <- X
  lnnxxvec <- log(X)
  ind <- nxxvec != 0
  nxxvec <- nxxvec[ind]
  lnnxxvec <- lnnxxvec[ind]
  LambdaA <- sum(nxxvec * lnnxxvec)
  nxvec <- c(nA, nB)
  lnnxvec <- log(nxvec)
  ind <- nxvec != 0
  nxvec <- nxvec[ind]
  lnnxvec <- lnnxvec[ind]
  LambdaB <- sum(nxvec * lnnxvec)
  lnLambda <- X[2] * log(2) + LambdaB - n * log(4 * n) - LambdaA
  Lambda <- exp(lnLambda)
  G2 <- -2 * lnLambda
  pval <- pchisq(G2, 1, lower.tail = FALSE)
  if (verbose) {
    cat("Likelihood ratio test for Hardy-Weinberg equilibrium\n")
    cat("G2 =", G2, "DF = 1", "p-value =", pval, "\n")
  }
} else { # x-linked marker
  if (length(X) != 5 | any(X < 0)) 
    stop("HWLratio: X is not a 5 by 1 non-negative count vector for an x-linked marker")
  if (any(!is.wholenumber(X))) {
     warning("Genotype counts are not integers, counts will be rounded.")
     X <- round(X, digits = 0)
  }

  X <- order.x(X)

  n <- sum(X)         
    
  nA <- X[1] + 2*X[3] + X[4]
  nB <- X[2] + 2*X[5] + X[4]
    
  nm <- X[1] + X[2]
  nf <- n - nm
  nt <- nA+nB

  nxxvec <- X # genotype counts
  lnnxxvec <- log(X)
  ind <- nxxvec != 0
  nxxvec <- nxxvec[ind]
  lnnxxvec <- lnnxxvec[ind]
  LambdaA <- -sum(nxxvec * lnnxxvec)
         
  nxvec <- c(nA, nB) # allele counts
  lnnxvec <- log(nxvec)
  ind <- nxvec != 0
  nxvec <- nxvec[ind]
  lnnxvec <- lnnxvec[ind]
  LambdaB <- sum(nxvec * lnnxvec)
         
  nsexvec <- c(nm,nf) # sex counts
  lnnsexvec <- log(nsexvec)
  ind <- nsexvec != 0
  nsexvec <- nsexvec[ind]    
  lnnsexvec <- lnnsexvec[ind]    
  LambdaC <- sum(nsexvec * lnnsexvec)
        
  lnLambda <- LambdaA + LambdaB + LambdaC + X[4]*log(2) - nt*log(nt)
  DF <- 2
    
  Lambda <- exp(lnLambda)
  G2 <- -2 * lnLambda
  pval <- pchisq(G2, DF, lower.tail = FALSE)
  if (verbose) {
    cat("Likelihood ratio test for Hardy-Weinberg equilibrium for an X-linked marker\n")
    cat("G2 =", G2, "DF =", DF, "p-value =", pval, "\n")
  }         
  } # else if x.linkd
out <- list(Lambda = Lambda, G2 = G2, pval = pval)
}
