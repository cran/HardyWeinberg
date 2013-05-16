HWLratio <- function (X, verbose = FALSE) 
{
    # Likelihood ratio test for HWE.
    if(length(X)!=3 | any(X<0)) stop("HWLratio: X is not a 3 by 1 non-negative count vector")
    if(any(!is.wholenumber(X))) {
      warning("Genotype counts are not integers, counts will be rounded.")
      X <- round(X,digits=0)
    }
    Xhom <- X[homozyg(X)]
    Xhet <- X[heterozyg(X)]
    X <- c(min(Xhom),Xhet,max(Xhom))  
    nAA <- X[1]
    nAB <- X[2]
    nBB <- X[3]
    n <- sum(X)
    nA <- 2 * nAA + nAB # allele counts
    nB <- 2 * nBB + nAB

    nxxvec <- X
    lnnxxvec <- log(X)
    ind <- nxxvec!=0
    nxxvec <- nxxvec[ind]
    lnnxxvec <- lnnxxvec[ind]
    LambdaA <- sum(nxxvec*lnnxxvec)

    nxvec <- c(nA,nB)
    lnnxvec <- log(nxvec)
    ind <- nxvec!=0
    nxvec <- nxvec[ind]
    lnnxvec <- lnnxvec[ind]
    LambdaB <- sum(nxvec*lnnxvec)
    
    lnLambda <- nAB * log(2) + LambdaB - n * log(4 * n) - LambdaA

    Lambda <- exp(lnLambda)
    G2 <- -2 * lnLambda
    pval <- pchisq(G2, 1, lower.tail = FALSE)
    if (verbose) {
        cat("Likelihood ratio test for Hardy-Weinberg equilibrium\n")
        cat("G2 =", G2, "p-value =", pval, "\n")
      }
    return(list(Lambda = Lambda, G2 = G2, pval = pval))
}

