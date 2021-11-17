HWChisq <- function (X, cc = 0.5, verbose = TRUE, x.linked = FALSE, phifixed = NULL) 
{
  if(!x.linked) { # autosomal
    if (is.matrix(X)) { #multiple alleles
      N <- sum(X)
      k <- ncol(X)
      ng <- 0.5*k*(k+1)
      ac <- rowSums(X)+colSums(X)
      alf <- ac/sum(ac)
      EP <- alf%o%alf
      dd <- diag(EP)
      EP <- 2*EP
      EP[upper.tri(EP)] <- 0
      diag(EP) <- dd
      He <- sum(EP[lower.tri(EP)])
      Ho <- sum(X[lower.tri(X)])/N
      EC <- EP*N
      D <- X-EC
      chi.contrib <- (D*D)/EC
      chi.contrib[upper.tri(chi.contrib)] <- 0
      chisq <- sum(chi.contrib)
      df <- 0.5*k*(k-1)
      pval <- pchisq(chisq,df=df,lower.tail=FALSE)
      D <- NA
      f <- NA
      p <- NA
      expected <- EC
      if (verbose) {
         title <- "Chi-square test for Hardy-Weinberg equilibrium (autosomal; multiple alleles)\n"
         cat(title)
	 cat(k, "alleles detected.\n")
         cat("Chi2 = ", chisq, "DF = ", df, "p-value = ", pval, "Ho = ", Ho, "He = ", He, "\n")
      }
    } else { #bi allelic
      if (length(X) != 3 | any(X < 0)) 
        stop("HWChisq: X is not a 3 by 1 non-negative count vector")
      if (any(!is.wholenumber(X))) {
        warning("Genotype counts are not integers, counts will be rounded.")
        X <- round(X, digits = 0)
      }
      mono <- FALSE
      DF <- 1
      ccv <- rep(cc, 3)
      n <- sum(X)
      Xhom <- X[homozyg(X)]
      Xhet <- X[heterozyg(X)]
      X <- c(min(Xhom), Xhet, max(Xhom))
      names(X) <- c("AA", "AB", "BB")
      p <- (2 * X[1] + X[2])/(2 * n)
      names(p) <- NULL
      q <- 1 - p
      if ((p == 0) | (q == 0)) {
        if (verbose) 
          cat("warning: monomorphic marker\n")
          mono <- TRUE
      }
      obs <- X
      expected <- c(n * p^2, n * 2 * p * q, n * q^2)
      names(expected) <- names(X)
      if (any(expected < 5)) 
        warning("Expected counts below 5: chi-square approximation may be incorrect")
      D <- 0.5 * (obs[2] - expected[2])
      names(D) <- NULL
      chi <- (abs(obs - expected) - ccv)^2
      f <- HWf(X)
      if (!mono) {
        chi.contrib <- chi/expected
        chisq <- sum(chi.contrib)
        pval <- pchisq(chisq, 1, lower.tail = FALSE)
      }
      else {
        chisq <- NA
	chi.contrib <- NA
        pval <- 1
      }
      if (verbose) {
        if (cc == 0) 
          title <- "Chi-square test for Hardy-Weinberg equilibrium (autosomal)\n"
        else title <- "Chi-square test with continuity correction for Hardy-Weinberg equilibrium (autosomal)\n"
        cat(title)
        cat("Chi2 = ", chisq, "DF = ", DF, "p-value = ", pval, "D = ", D, 
            "f = ", f, "\n")
      }
    } # bi-allelic
  } else { # sex-linked
      if (length(X) != 5 | any(X < 0)) 
        stop("HWChisq: X is not a 5 by 1 non-negative count vector for a X-linked marker")
      if (any(!is.wholenumber(X))) {
        warning("Genotype counts are not integers, counts will be rounded.")
        X <- round(X, digits = 0)
      }
      mono <- FALSE
      ccv <- rep(cc, 5)
      n <- sum(X)
      lab <- names(X)
      if(!all(lab %in% c("A","AA","AB","B","BB")))
        stop("Unknown genotypes occurred. Supply counts as a named vector c(A,AA,AB,B,BB)")
      nfAA <- X[lab=="AA"]
      nfAB <- X[lab=="AB"]
      nfBB <- X[lab=="BB"]
      nmA <- X[lab=="A"]
      nmB <- X[lab=="B"]
      #    mcat(nfAA,nfAB,nfBB)
      #    mcat(nmA,nmB)
      nm <- nmA+nmB
      nf <- n - nm
      X <- c(nmA,nmB,nfAA,nfAB,nfBB)
      names(X) <- c("A","B","AA", "AB", "BB")
      p <- (2 * nfAA + nfAB + nmA)/(2 * nf + nm)
      if(is.null(phifixed)) theta <- nm/n else theta <- phifixed
      #    mcat(theta,p)
      names(p) <- NULL
      q <- 1 - p
      if ((p == 0) | (q == 0)) {
        if (verbose) 
          cat("warning: monomorphic marker\n")
        mono <- TRUE
      }
      obs <- X
      expected <- c(theta*n*p, theta*n*(1-p), (1-theta)*n * p^2, (1-theta)*n * 2 * p * q, (1-theta)*n * q^2)
      names(expected) <- names(X)
      if (any(expected < 5)) 
        warning("Expected counts below 5: chi-square approximation may be incorrect")
      chi <- (abs(obs - expected) - ccv)^2
      if(all(X[3:5]==0)) f <- NA else f <- HWf(X[3:5]) # computed from females only.
      D <- NA
      if(is.null(phifixed)) DF <- 2 else DF <- 3
      if (!mono) {
	chi.contrib <- chi/expected
        chisq <- sum(chi.contrib)    
        pval <- pchisq(chisq, DF, lower.tail = FALSE)
      }
      else {
        chisq <- NA
	chi.contrib <- NA
        pval <- 1
      }
      if (verbose) {
        if (cc == 0) 
          title <- "Chi-square test for Hardy-Weinberg equilibrium (X-chromosomal)\n"
        else title <- "Chi-square test with continuity correction for Hardy-Weinberg equilibrium (X-chromosomal)\n"
        cat(title)
        cat("Chi2 = ", chisq, "DF =", DF, "p-value = ", pval, "D = ", D, 
            "f = ", f, "\n")
      }
  } # end else !x.linked
  out <- list(chisq = chisq, pval = pval, D = D, p = p, f = f, expected = expected, chi.contrib = chi.contrib)
}
