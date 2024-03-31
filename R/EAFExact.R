EAFExact <- function(m,f,verbose=TRUE,...) {
    knownformat <- FALSE
    m <- unlist(m)
    f <- unlist(f)
  nm <- sum(m)
  nf <- sum(f)
  if(is.vector(m)) {
      alleles.m <- alleles(m)
  } else {
      alleles.m <- colnames(m)
  }
  if(is.vector(f)) {
    alleles.f <- alleles(f)
  } else {
    alleles.f <- colnames(f)
  }
  allele.names.equal <- all(sort(alleles.m) == sort(alleles.f))
  if(!allele.names.equal) stop("EAF: male and female alleles different")
  if(is.vector(m) & is.vector(f)) {
    if(length(m)==2 & length(f)==3) {
      # X-chromosomal bi-alleliec
      knownformat <- TRUE
      lab.f <- names(f)
      hom1 <- paste(alleles.f[1],alleles.f[1],sep="")
      hom2 <- paste(alleles.f[2],alleles.f[2],sep="")
      het  <- names(f)[heterozyg(f)]
      nfAA <- f[lab.f == hom1]
      nfAB <- f[lab.f == het]
      nfBB <- f[lab.f == hom2]
      fA <- 2 * nfAA + nfAB
      fB <- 2 * nfBB + nfAB
      lab.m <- names(m)
      nmA <- m[lab.m == alleles.f[1]]
      nmB <- m[lab.m == alleles.f[2]]
      m.ac <- c(nmA, nmB)
      f.ac <- c(fA, fB)
      nt <- 2 * nf + nm
    }
    if(length(m)==3 & length(f)==3) {
      # Autosomomal bi-allelic
      knownformat <- TRUE
      lab.f <- names(f)
      hom1 <- paste(alleles.f[1],alleles.f[1],sep="")
      hom2 <- paste(alleles.f[2],alleles.f[2],sep="")
      hetf  <- names(f)[heterozyg(f)]
      hetm  <- names(m)[heterozyg(m)]
      if(hetf!=hetm) {
        stop("heterozygotes of males and females not equally coded.\n")
      }
      nfAA <- f[lab.f == hom1]
      nfAB <- f[lab.f == hetf]
      nfBB <- f[lab.f == hom2]
      fA <- 2 * nfAA + nfAB
      fB <- 2 * nfBB + nfAB
      lab.m <- names(m)
      nmAA <- m[lab.m == hom1]
      nmAB <- m[lab.m == hetf]
      nmBB <- m[lab.m == hom2]
      mA <- 2 * nmAA + nmAB
      mB <- 2 * nmBB + nmAB
      m.ac <- c(mA, mB)
      f.ac <- c(fA, fB)
      nt <- 2 * (nm + nf)
    }
    if(length(m) > 3 & length(f) > 3) { # these will be counted as autosomal multiallelic below.
      m <- toTriangular(m)
      f <- toTriangular(f)
    }
  } 
  if(is.vector(m) & is.matrix(f)) {
    # X chromosomal multi allelic
    knownformat <- TRUE
    m.ac <- m
    f.ac <- rowSums(f) + colSums(f)
    nt <- sum(m) + 2*sum(f)
  }
  if(is.matrix(m) & is.matrix(f)) {
    # autosomal multiallelic
    knownformat <- TRUE
    f.ac <- rowSums(f) + colSums(f)
    m.ac <- rowSums(m) + colSums(m)
    nt <- 2*(sum(m) + sum(f))
  }
  if(!knownformat) stop("unknown format, revise input of m and f")
  tab <- rbind(m.ac,f.ac)
  rownames(tab) <- c("males","females")
  colnames(tab) <- alleles.f
  out <- fisher.test(tab,...)
  pval <- out$p.value
  n <- nm + nf
  if(verbose) {
    cat("Fisher Exact test for equality of allele frequencies for males and females.\n\n")
    cat("Table of allele counts:\n\n")
    print(tab)
    cat("\nSample of",n,"indivduals with ",nt,"alleles. p.value = ",pval,"\n")
  }
  out <- list(pval=pval,tab=tab)
}
