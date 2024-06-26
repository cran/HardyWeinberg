HWExact <- function (X, alternative = "two.sided", pvaluetype = "selome", eps=1e-10, x.linked = FALSE, verbose = TRUE)
{
  if(!x.linked) {
    if (length(X) != 3 | any(X < 0))
      stop("HWExact: X is not a 3 by 1 non-negative count vector")
    if (any(!is.wholenumber(X))) {
      warning("Genotype counts are not integers, counts will be rounded.")
      X <- round(X, digits = 0)
    }
    n <- sum(X)
    Xhom <- X[homozyg(X)]
    Xhet <- X[heterozyg(X)]
    nA <- 2 * Xhom[1] + Xhet
    nB <- 2 * n - nA
    MaxHet <- min(nA, nB)
    if (MaxHet < 2) {
      pval <- 1
      prob <- 1
      pofthesample <- 1
      ind <- 1
    }
    else {
      ind <- match(Xhet, seq(MaxHet%%2, MaxHet, 2))
      enAB <- nA * nB/(2 * n - 1)
      enAB <- round(enAB, digits = 0)
      if ((enAB%%2) != (MaxHet%%2))
        enAB <- enAB + 1
      nAA <- (nA - enAB)/2
      nBB <- (nB - enAB)/2
      initialprob <- 1
      AboveExp <- NULL
      BelowExp <- NULL
      if (enAB < MaxHet)
        AboveExp <- CompProbUp(nAA, nBB, enAB, initialprob,
                               MaxHet)
      BelowExp <- CompProbDown(nAA, nBB, enAB, initialprob)
      prob <- c(rev(BelowExp), initialprob, AboveExp)
      prob <- prob/sum(prob)
    }
    if (MaxHet%%2 == 0)
      names(prob) <- seq(0, MaxHet, 2)
    if (MaxHet%%2 == 1)
      names(prob) <- seq(1, MaxHet, 2)
    Plow <- cumsum(prob)
    Phigh <- 1 - c(0, Plow)
    Phigh <- Phigh[-length(Phigh)]
    Phwe <- pmin(1, 2 * Phigh, 2 * Plow)
    pofthesample <- prob[ind]
    pval <- switch(alternative,
                   greater = switch(pvaluetype, selome = Phigh[ind], midp = Phigh[ind] - 0.5 * pofthesample,
                                    stop("invalid value for parameter pvaluetype")),
                   less = switch(pvaluetype, selome = Plow[ind], midp = Plow[ind] - 0.5 * pofthesample,
                                 stop("invalid value for parameter pvaluetype")),
                   two.sided = switch(pvaluetype, dost = Phwe[ind], selome = sum(prob[prob <=
                                                                                        pofthesample]), midp = sum(prob[prob < pofthesample]) +
                                        0.5 * pofthesample, stop("invalid value for parameter pvaluetype")),
                   stop("invalid value for parameter alternative"))
    if (verbose) {
      D <- 0.5 * (Xhet - nA * nB/(2 * n))
      cat("Haldane Exact test for Hardy-Weinberg equilibrium (autosomal)\n")
      stringpvalue <- switch(pvaluetype, dost = "using DOST p-value\n",
                             selome = "using standard (SELOME) p-value\n", midp = "using MID p-value\n",
                             stop("invalid value for parameter pvaluetype"))
      cat(stringpvalue)
      cat(paste("sample counts: n", names(Xhom[1]), " = ",
                sep = ""), Xhom[1], paste("n", names(Xhet), " = ",
                                          sep = ""), Xhet, paste("n", names(Xhom[2]), " = ",
                                                                 sep = ""), Xhom[2], "\n")
      stringtwosided <- paste("H0: HWE (D==0), H1: D <> 0 \nD = ",
                              format(D, scientific = FALSE), "p-value = ", format(pval,
                                                                            scientific = FALSE), "\n")
      stringgreater <- paste("H0: HWE (D==0), H1: D > 0 \nD = ",
                             format(D, scientific = FALSE), "p-value = ", format(pval,
                                                                           scientific = FALSE), "\n")
      stringless <- paste("H0: HWE (D==0), H1: D < 0 \nD = ",
                          format(D, scientific = FALSE), "p-value = ", format(pval,
                                                                        scientific = FALSE), "\n")
      toprint <- switch(alternative, two.sided = stringtwosided,
                        greater = stringgreater, less = stringless)
      cat(toprint)
    }
  } else { # x.linked marker.
    if (length(X) != 5 | any(X < 0)) 
      stop("HWExact: X is not a 5 by 1 non-negative count vector for an x-linked marker")
    if (any(!is.wholenumber(X))) {
      warning("Genotype counts are not integers, counts will be rounded.")
      X <- round(X, digits = 0)
    }
    lab <- names(X)
    if(is.null(lab)) {
      warning("Genotype counts not labelled, default sequence c(A,B,AA,AB,BB) assumed.")
      names(X) <- c("A","B","AA","AB","BB")
    }
    als <- alleles(X)
    if(length(als) != 2) {
       stop("HWExact: X is not a bi-allelic x-linked marker")
    }
    hom1 <- paste(als[1],als[1],sep="")
    hom2 <- paste(als[2],als[2],sep="")
    fems <- X[!(names(X) %in% als)]
    het  <- names(fems)[heterozyg(fems)]
    theo.types <- c(als[1],als[2],hom1,het,hom2)
    if(!all(lab %in%  theo.types))
      stop("Unknown genotypes occurred. Supply counts as a named vector c(A,B,AA,AB,BB)")
    n <- sum(X)         

    X <- order.x(X)
   
    ### recompute all genotype and allele counts considering A the minor allele

    nfAA <- X[3]
    nfAB <- X[4]
    nfBB <- X[5]
    nmA <- X[1]
    nmB <- X[2]

    nm <- nmA+nmB
    nf <- n - nm

    nAf <- 2*nfAA + nfAB
    nBf <- 2*nfBB + nfAB
    nA <- nmA + 2*nfAA + nfAB
    nB <- nmB + 2*nfBB + nfAB
    
    nt <- nA+nB
    pA <- nA/nt
   
    Z <- auxiliartable(X)
    
    prob <- numeric(nrow(Z))
    
    pofthesample <- sample.prob.last(n,nm,nmA,nA,nfAB)

    for(i in 1:nrow(Z)) {
      prob[i] <- subsamples.prob(nA, nB, nm, nf, nt, Z[i,1], Z[i,2], Z[i,3], pofthesample, eps)
    }

    if(pvaluetype=="selome") pval <- sum(prob)
    if(pvaluetype=="midp") {
      pval <- sum(prob)-0.5*pofthesample
    }
    if(verbose) {
      cat("Graffelman-Weir exact test for Hardy-Weinberg equilibrium on the X-chromosome\n")
      stringpvalue <- switch(pvaluetype, 
                             selome = "using SELOME p-value\n", midp = "using MID p-value\n",
                             stop("invalid value for parameter pvaluetype"))
      cat(stringpvalue)
      cat("Sample probability",pofthesample,"p-value = ",pval,"\n")
    }
    prob <- NA # there is no final probability vector in the x.linked case
  }
  out <- list(pval = pval, prob = prob, pofthesample = pofthesample)
}
