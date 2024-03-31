HWStr <- function (X, test = "permutation", verbose = FALSE, ...) 
{
    nc <- ncol(X)
    n  <- nrow(X)
    nstrs <- nc/2
    cat(nstrs, "STRs detected.\n")
    a1 <- seq(1, nc, 2)
    a2 <- seq(2, nc, 2)
    STR <- colnames(X)[a1]
    nr <- 1:nstrs
    Results <- data.frame(STR)
    ss    <- numeric(nstrs)
    pval  <- numeric(nstrs)
    minor <- numeric(nstrs)
    major <- numeric(nstrs)
    nals  <- numeric(nstrs)
    ho    <- numeric(nstrs)
    he    <- numeric(nstrs)
    Hp    <- numeric(nstrs)
    for (i in 1:nstrs) {
      imis <- is.na(X[, a1[i]]) | is.na(X[, a1[i]])
      nmis <- sum(imis)
      ss[i] <- n - nmis
      GC <- AllelesToTriangular(X[, a1[i]], X[, a2[i]])
      if(verbose) print(GC)
      nals[i] <- nrow(GC)
      ac <- colSums(GC) + rowSums(GC)
      nt <- sum(ac)
      twon <- 2*ss[i]
      if(nt!=twon) {
        stop("inconsistent allele counts.\n")
      }
      af  <- ac/nt
      saf <- sum(af)
      if (saf!=1) {
        stop("allele frequencies do not sum 1.\n")
      }
      minor[i] <- min(ac)/nt
      major[i] <- max(ac)/nt
      ho[i] <- (n - sum(diag(GC)))/n  # observed heterozygosity
      he[i] <- 1-sum(af^2)            # expected heterozygosity
      Hp[i] <- shannon(ac)$Hp         # shannon index of allele frequencies
      if (test == "none") {
          pval[i] <- NA
      } else if (test == "chisq") {
        out <- HWChisq(GC, cc = 0, verbose = verbose)
        pval[i] <- out$pval
      }
      else if (test == "permutation") {
        out <- HWPerm.mult(GC, verbose = FALSE, ...)
        pval[i] <- out$pval
      }
      else {
        cat("Unknown type of test. P-values set to NA.\n")
        pval[i] <- NA
      }
    }
    maxminoraf <- 1/nals
    check <- minor > maxminoraf
    if(any(check)) {
      cat("minor allele frequency out of range\n")
      for(j in 1:nstrs) {
        if(check[j]) {
          cat(j,maxminoraf[j],minor[j],"\n")    
        }
      }
    }
    maxhe <- (nals-1)/nals
    check <- he > maxhe
    if(any(check)) {
      cat("expected heterozygosity out of range\n")
      for(j in 1:nstrs) {
        if(check[j]) {
          cat(j,maxhe[j],he[j],"\n")    
        }
      }
    }
    Results$N <- ss
    Results$Nt <- nals
    Results$MinorAF <- minor
    Results$MajorAF <- major
    Results$Ho <- ho
    Results$He <- he
    Results$Hp <- Hp
    Results$pval <- pval
    return(Results)
}
