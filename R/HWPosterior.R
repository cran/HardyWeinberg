HWPosterior <- function(males, females, verbose = TRUE, prior.af = c(0.5, 0.5), 
                            prior.gf = c(0.333, 0.333, 0.333), x.linked = FALSE, 
                            precision = 0.05) 
{
  if(x.linked) { # X-chromosomal
    X <- c(males, females)
    if (length(X) != 5) {
      stop("Improper number of genotype counts.")
    }
    X <- order.x(X)
    n <- sum(X)
    x.m <- X[1:2]
    x.f <- X[3:5]
    p.h0 <- H0(x.f, x.m, alpha = prior.af)
    p.h1 <- H1(x.f, x.m, alpha = prior.gf)
    p.h2 <- H2(x.f, x.m, alpha.f = prior.af, alpha.m = prior.af)
    p.h3 <- H3(x.f, x.m, alpha.f = prior.gf, alpha.m = prior.af)
    posterior.prob <- c(p.h0, p.h1, p.h2, p.h3)
    posterior.prob <- posterior.prob/sum(posterior.prob)
    lBF <- log10(3 * posterior.prob/(1 - posterior.prob))
    r.posterior.prob <- round(posterior.prob, 4)
    Res <- cbind(posterior.prob, lBF)
    colnames(Res) <- c("Posterior_Prob", "log10(Bayes Factor)")
    rownames(Res) <- c("M0 (HWE):", "M1 (f!=0):", "M2 (d!=1):", 
                       "M3 (f!=0 & d!=1:)")
    if (verbose) {
      cat("Bayesian test for Hardy-Weinberg equilibrium of X-chromosomal variants.\n\n")
      print(round(Res, digits = 4))
    }
    out <- Res 
  } else { # autosomal

    if (length(males) != 3 | length(females) != 3) {
      stop("Improper numbers of genotype counts.")
    }
    
    males   <- order.auto(males)
    females <- females[names(males)] # line up the counts.
    total   <- males + females
    total   <- order.auto(total) # sequence according to overall minor allele
    males   <- males[names(total)]
    females <- females[names(total)]

    X <- c(males,females)

    posterior.prob <- c(M11p(X,alpha=prior.af),
              M12p(X,alpha=prior.gf,r=precision),
              M13p(X,alpha=prior.gf,r=precision),
              M14p(X,alpha=prior.gf),
              M15p(X,alpha.f=prior.gf,r=precision), # alpha.m default
              M21p(X,alpha.f=prior.af,alpha.m=prior.af),
              M22p(X,alpha.f=prior.gf,alpha.m=prior.af),
              M23p(X,alpha.f=prior.af,alpha.m=prior.gf),
              M24p(X,alpha.f=prior.gf,alpha.m=prior.af,r=precision),
              M25p(X,alpha.f=prior.gf,alpha.m=prior.gf))
    posterior.prob <- posterior.prob/sum(posterior.prob)
    names(posterior.prob) <- c("M_11", "M_12", "M_13", "M_14", "M_15", "M_21", 
                    "M_22", "M_23", "M_24", "M_25")
    if(verbose) {
      print(round(posterior.prob,digits=4))
      indexbest <- which.max(posterior.prob)
      cat("Best fitting", names(posterior.prob[indexbest]), posterior.prob[indexbest], 
          "\n")
    }
    out <- posterior.prob
  }
}
