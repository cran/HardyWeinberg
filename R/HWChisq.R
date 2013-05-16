HWChisq <- function (X, cc = 0.5, alpha = 0.05, verbose = FALSE) 
{
    if(length(X)!=3 | any(X<0)) stop("HWChisq: X is not a 3 by 1 non-negative count vector")
    if(any(!is.wholenumber(X))) {
      warning("Genotype counts are not integers, counts will be rounded.")
      X <- round(X,digits=0)
    }
    if(alpha <= 0 | alpha >=1) stop("HWChisq: alpha should be in the range (0,1)")
    mono <- FALSE
    ccv <- rep(cc, 3)
    n <- sum(X)
    Xhom <- X[homozyg(X)]
    Xhet <- X[heterozyg(X)]
    # reorganize counts in standard AA, AB, BB format.
    X <- c(min(Xhom),Xhet,max(Xhom))
    names(X) <- c("AA","AB","BB")
    p <- (2 * X[1] + X[2])/(2 * n)
    names(p) <- NULL
    q <- 1 - p
    if ((p == 0) | (q == 0)) {
        if (verbose) 
            cat("warning: monomorphic marker\n")
        mono <- TRUE
    }
    obs <- X
    # calculate expected counts under HWE.
    exp <- c(n * p^2, n * 2 * p * q, n * q^2)
    D <- 0.5 * (obs[2] - exp[2])
    names(D) <- NULL
    # determine a cutoff value for the minor allele frequency, below
    # this cutoff, a continuity correction is no longer applied.
#    root <- cutoff(n,alpha,cc,verbose=FALSE)
#    if(maf(X)<root) {
#      cc <- 0
#      ccv <- rep(0, 3)
#    }
    chi <- (abs(obs - exp) - ccv)^2
    # compute the inbreeding coefficient
    f <- HWf(X)
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
       cat("Chi2 = ", chisq, "p-value = ", pval, "D = ", D,"f = ", f, 
            "\n")
    }
    return(list(chisq = chisq, pval = pval, D = D, p = p, f = f))
}

