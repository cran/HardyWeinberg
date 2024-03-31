HWData <- function (nm = 100, n = rep(100, nm), f = rep(0, nm), p = NULL, 
    conditional = FALSE, exactequilibrium = FALSE, x.linked = FALSE, nA = NULL, 
    n.males = round(0.5 * n), shape1 = 1, shape2 = 1, counts = TRUE) {
  
  # total numbers of alleles
  
  n.females <- n - n.males
  
  if(!x.linked) {
    nt <- 2*n
  } else {
    nt <- n.males + 2*n.females
  }
  
  # uniform allele frequencies if unspecified
  
  if(is.null(p) & is.null(nA)) { # nothing specified
    p <- rbeta(nm, shape1 = shape1, shape2 = shape2)
  } else if(is.null(p)) { # only counts given
    p <- nA/nt
  } else { # only p given
    nA <- as.integer(round(p*nt,digits=0))
  }
  
  # extend scalars to vectors
  
  if (length(p) == 1) p <- rep(p, nm)
  if (length(n) == 1) n <- rep(n, nm)
  if (length(f) == 1) f <- rep(f, nm)
  if (length(nA) == 1) nA <- rep(nA, nm)
  if (length(n.males) == 1) n.males <- rep(n.males, nm)
  
  # output matrix
  
  if(!x.linked) {
    X <- matrix(0, nrow = nm, ncol = 3)
    colnames(X) <- c("AA", "AB", "BB")   
  } else {
    X <- matrix(0, nrow = nm, ncol = 5)
    colnames(X) <- c("A","B","AA", "AB", "BB")   
  }
  
  # generate genotype counts
  
  if(!x.linked) { # autosomal
    if(conditional) { # levene-haldane
      if(exactequilibrium) {
        nAA <- round(n*p*p,digits=0)
        nAB <- round(2*n*p*(1-p),digits=0)
        nBB <- round(n*((1-p)^2),digits=0)
        X <- cbind(nAA,nAB,nBB)
	
      } else {
        nB <- nt - nA
        for (i in 1:nm) {
          Pop <- c(rep(1, nA[i]), rep(0, nB[i]))
          sam <- matrix(sample(Pop), ncol = 2, byrow = TRUE)
          status <- apply(sam, 1, sum)
          nAA <- sum(status == 2)
          nAB <- sum(status == 1)
          nBB <- sum(status == 0)
          if ((nAA + nAB + nBB) != n[i]) stop("HWData: error")
          X[i, ] <- c(nAA, nAB, nBB)
        }
      }
    } else { # multinomial
      if(exactequilibrium) {
        X <- cbind(n*p^2 + n*f*p*(1 - p), n*(1 - f)*2*p*(1 - p), 
              n*(1 - p)^2 + n*f*p*(1 - p) ) 
        if(counts) {
          X <- round(X,digits=0)
        }
      } else {
        for (i in 1:nm) {
          X[i, ] <- t(rmultinom(1, size = n[i], prob = c(p[i]^2 + 
                      f[i] * p[i] * (1 - p[i]), (1 - f[i]) * 2 * 
                      p[i] * (1 - p[i]), (1 - p[i])^2 + f[i] * 
                      p[i] * (1 - p[i]))))
        }  
      } # end if exact equilibrium
    }
    colnames(X) <- c("AA", "AB", "BB")
  } else { # X-chromosomal
    theta <- n.males/n
    if(conditional) { #graffelman-weir
      for (i in 1:nm) {
        Pop <- c(rep(1, nA[i]), rep(0, nt[i] - nA[i]))
        Sam <- sample(Pop)
        Males <- Sam[1:n.males[i]]
        mA <- sum(Males == 1)
        mB <- sum(Males == 0)
        Females <- matrix(Sam[(n.males[i] + 1):nt[i]], ncol = 2, byrow = TRUE)
        status <- apply(Females, 1, sum)
        fAA <- sum(status == 2)
        fAB <- sum(status == 1)
        fBB <- sum(status == 0)
        if ((fAA + fAB + fBB) != n.females[i]) stop("HWData: error")
        X[i, ] <- c(mA, mB, fAA, fAB, fBB)
      }
    } else { #multinomial
      if(exactequilibrium) {
        X <- cbind(n.males*theta*p, 
                   n.males*theta*(1 - p), 
                   n.females*(1 - theta) * p^2, 
                   n.females*(1 - theta) * 2 * p * (1 - p), 
                   n.females*(1 - theta) * (1 - p)^2)
        if(counts) {
          X <- round(X,digits=0)
        }
      } else {
        prob <- cbind(theta*p, 
                    theta*(1 - p), 
                    (1 - theta) * p^2, 
                    (1 - theta) * 2 * p * (1 - p), 
                    (1 - theta) * (1 - p)^2)
        for (i in 1:nm) {
          X[i,] <- t(rmultinom(1, size = n[i], prob[i, ]))
        }
      } # end if exactequilibrium
    }
    colnames(X) <- c("A","B","AA", "AB", "BB")   
  } # end if x-chromosomal
  if(!counts) {
    X <- X/rowSums(X)
  }
  return(X)
} 

