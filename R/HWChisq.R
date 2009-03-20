HWChisq <-
function (X, cc = 0.5, alpha = 0.05, verbose = FALSE) 
{
    mono <- FALSE
    ccv <- rep(cc, 3)
    n <- sum(X)
    p <- (2 * X[1] + X[2])/(2 * n)
    q <- 1 - p
    if ((p == 0) | (q == 0)) {
        if (verbose) 
            cat("warning: monomorphic marker\n")
        mono <- TRUE
    }
    obs <- X
    exp <- c(n * p^2, n * 2 * p * q, n * q^2)
    D <- 0.5 * (obs[2] - exp[2])
    root <- cutoff(n,alpha,cc,verbose=FALSE)
    if(maf(X)<root) {
      cc <- 0
      ccv <- rep(0, 3)
    }
    chi <- (abs(obs - exp) - ccv)^2
    if (!mono) {
        chi2 <- chi/exp
        chisq <- sum(chi2)
        pval <- 1 - pchisq(chisq, 1)
    }
    else {
        chisq <- NA
        pval <- 1
    }
    if (verbose) {
       if(cc==0) title <- "Chi-square test for Hardy-Weinberg equilibrium\n"
       else title <- "Chi-square test with continuity correction for Hardy-Weinberg equilibrium\n"
       cat(title)
       cat("Chi2 = ", chisq, "p-value = ", pval, "D = ", D, 
            "\n")
    }
    return(list(chisq = chisq, pval = pval, D = D, p = p))
}

