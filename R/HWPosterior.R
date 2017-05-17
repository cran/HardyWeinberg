HWPosterior <- function(X,verbose=TRUE,prior.af=c(0.5,0.5),prior.gf=c(0.333,0.333,0.333)) {
  lab <- names(X)
  if(is.null(lab)) {
    cat("No genotype labels given, default order c(A,B,AA,AB,BB) assumed.\n")
    X <- genlabels(X)
    lab <- names(X)
  }
  if(length(X)!=5) {
    cat("Improper number of genotype counts.\n")
    stop()
  }
  if (!all(lab %in% c("A", "B", "AA", "AB", "BB"))) 
    stop("Unknown genotypes occurred. Supply counts as a named vector c(A,AA,AB,B,BB)")
  n <- sum(X)
  nfAA <- X[lab == "AA"]
  nfAB <- X[lab == "AB"]
  nfBB <- X[lab == "BB"]
  nmA <- X[lab == "A"]
  nmB <- X[lab == "B"]
  x.m <- c(nmA,nmB)
  x.f <- c(nfAA,nfAB,nfBB)
  p.h0 <- H0(x.f, x.m, alpha = prior.af)
  p.h1 <- H1(x.f, x.m, alpha = prior.gf) 
  p.h2 <- H2(x.f, x.m, alpha.f = prior.af, alpha.m = prior.af)
  p.h3<- H3(x.f, x.m, alpha.f = prior.gf, alpha.m = prior.af)
  posterior.prob <- c(p.h0,p.h1,p.h2,p.h3)
  posterior.prob <- posterior.prob/sum(posterior.prob)
  lBF <- log10(3*posterior.prob/(1-posterior.prob))
  r.posterior.prob <- round(posterior.prob,4)
  Res <- cbind(posterior.prob,lBF)
  colnames(Res) <- c("Posterior_Prob","log10(Bayes Factor)")
  rownames(Res) <- c("M0 (HWE):","M1 (f!=0):",
                     "M2 (d!=1):","M3 (f!=0 & d!=1:)")
  if(verbose) {
    cat("Bayesian test for Hardy-Weinberg equilibrium of X-chromosomal variants.\n\n")
    print(round(Res,digits=4))
  }
  out <- Res
}
