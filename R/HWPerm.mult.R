HWPerm.mult <- function (x, y = NULL, nperm = 17000, eps = 1e-10, verbose = TRUE, ...) 
{
  testtype <- NULL
  if(is.vector(x) & is.null(y)) {
    x <- toTriangular(x)
    testtype <- 1 # autosomal
  } else if(is.matrix(x) & is.null(y)) {
    testtype <- 1 # autosomal
  } else if(is.matrix(x) & is.matrix(y)) {
    testtype <- 2 # autosomal, stratified by sex
  } else if(is.vector(x) & is.matrix(y)) {
    testtype <- 3 # X-chromosomal
  } else {
    stop("data vectors x and/or y not correctly specified.")
  }
  if (testtype == 1) { # autosomal
    n <- sum(x)
    k <- nrow(x)
    if (k != ncol(x)) 
      stop("error x is not square.")
    pofthesample <- dlevene(x)
    als <- colnames(x)
    ac.counts <- rowSums(x) + colSums(x)
    alvec <- rep(als, times = ac.counts)
    pseudodist <- numeric(nperm)
    for (l in 1:nperm) {
      scrambled <- sample(alvec)
      ind1 <- seq(1, 2 * n, 2)
      ind2 <- seq(2, 2 * n, 2)
      A1 <- scrambled[ind1]
      A2 <- scrambled[ind2]
      gmatrix.n <- AllelesToTriangular(A1,A2)
      pseudodist[l] <- dlevene(gmatrix.n)
    }
    ii <- nearlyEqual(pseudodist, rep(pofthesample, nperm), 
                                      eps)
    nsmaller <- sum(ii)
    iii <- !ii & pseudodist < pofthesample
    nsmaller <- nsmaller + sum(iii)
    pval <- nsmaller/nperm
    if (verbose) {
      cat("Permutation test for Hardy-Weinberg equilibrium (autosomal).\n")
      cat(k, "alleles detected.\n")
      cat("Observed statistic:", pofthesample, " ", 
          nperm, "permutations. p-value:", pval, 
          "\n")
    }
  }
  if (testtype == 2) { # autosomal, stratified by sex
    nm <- sum(x)
    nf <- sum(y)
    k <- nrow(x)
    m <- x
    f <- y
    pofthesample <- dens.auto(m, f)
    alsnames.m <- colnames(m)
    alsnames.f <- colnames(f)
    if (any(alsnames.m != alsnames.f)) 
      stop("different male and female alleles")
    mm.counts <- rowSums(m) + colSums(m)
    als.m <- rep(alsnames.m, times = mm.counts)
    ff.counts <- rowSums(f) + colSums(f)
    als.f <- rep(alsnames.f, times = ff.counts)
    alvec <- c(als.m, als.f)
    knownals <- sort(unique(alvec))
    knownals <- paste("A",knownals,sep="")
    pseudodist <- numeric(nperm)
    for (l in 1:nperm) {
      alvec.scrambled <- sample(alvec)
      males <- alvec.scrambled[1:(2 * nm)]
      females <- alvec.scrambled[(2 * nm + 1):length(alvec)]
      ind1.m <- seq(1, 2 * nm, 2)
      ind2.m <- seq(2, 2 * nm, 2)
      A1m <- males[ind1.m]
      A2m <- males[ind2.m]
      gmatrix.m <- AllelesToTriangular(A1m,A2m,given=knownals)
   
      ind1.f <- seq(1, 2 * nf, 2)
      ind2.f <- seq(2, 2 * nf, 2)
      A1f <- females[ind1.f]
      A2f <- females[ind2.f]
      gmatrix.f <- AllelesToTriangular(A1f,A2f,given=knownals)
      
      pseudodist[l] <- dens.auto(gmatrix.m, gmatrix.f)
    }
    ii <- nearlyEqual(pseudodist, rep(pofthesample, nperm), 
                      eps)
    nsmaller <- sum(ii)
    iii <- !ii & pseudodist < pofthesample
    nsmaller <- nsmaller + sum(iii)
    pval <- nsmaller/nperm
    if (verbose) {
      cat("Permutation test for Hardy-Weinberg equilibrium and equality of allele frequencies (autosomal).\n")
      cat(k, "alleles detected.\n")
      cat("Observed statistic:", pofthesample, " ", 
          nperm, "permutations. p-value:", pval, 
          "\n")
    }
  }
  if (testtype == 3) { # X-chromosomal
    nm <- sum(x)
    nf <- sum(y)
    m <- x
    f <- y
    k <- length(m)
    pofthesample <- density.ma.gender(m, f)
    als <- names(m)
    als.m <- rep(als, times = m)
    ff.counts <- rowSums(f) + colSums(f)
    als.f <- rep(als, times = ff.counts)
    alvec <- c(als.m, als.f)
    knownals <- sort(unique(alvec))
    knownals <- paste("A",knownals,sep="")
    pseudodist <- numeric(nperm)
    for (l in 1:nperm) {
      alvec.scrambled <- sample(alvec)
      males <- alvec.scrambled[1:nm]
      females <- alvec.scrambled[(nm + 1):length(alvec)]
      ind1.f <- seq(1, 2 * nf, 2)
      ind2.f <- seq(2, 2 * nf, 2)
      
      A1f <- females[ind1.f]
      A2f <- females[ind2.f]
      gmatrix.f <- AllelesToTriangular(A1f,A2f,given=knownals)
#      print(gmatrix.f)
      
      m.counts <- numeric(k)
      for (i in 1:length(als)) {
        m.counts[i] <- sum(males == als[i])
      }
      
      pseudodist[l] <- density.ma.gender(m.counts, gmatrix.f)
    }
    ii <- nearlyEqual(pseudodist, rep(pofthesample, nperm), 
                      eps)
    nsmaller <- sum(ii)
    iii <- !ii & pseudodist < pofthesample
    nsmaller <- nsmaller + sum(iii)
    pval <- nsmaller/nperm
    if (verbose) {
      cat("Permutation test for Hardy-Weinberg equilibrium and equality of allele frequencies (X-chromosomal).\n")
      cat(k, "alleles detected.\n")
      cat("Observed statistic:", pofthesample, " ", 
          nperm, "permutations. p-value:", pval, 
          "\n")
    }
  }
  out <- list(pofthesample = pofthesample, pseudodist = pseudodist, 
              pval = pval)
}
