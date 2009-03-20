HWData <-
function (n = 100, nm = 100, f = 0, p = NULL, pfixed = FALSE, exactequilibrium = FALSE, 
    pdist = "runif", ...) 
{
    Xt <- NULL
    if(!pfixed) { # multinomial sampling
      for (i in 1:nm) { 
         if (is.null(p)) {
            if (pdist == "runif") 
               ps <- runif(1, ...)
            else {
               if (pdist == "rbeta") 
                 ps <- rbeta(1, ...)
               else stop("unknown value for pdist")
            }
         }
         else
           ps <- p
         if (!exactequilibrium) 
           X <- t(rmultinom(1, size = n, prob = c(ps^2 + f*ps*(1-ps), (1-f)*2 * ps * 
               (1 - ps), (1 - ps)^2 + f*ps*(1-ps))))
         else X <- c(n * ps^2 + n*f*ps*(1-ps), n * (1-f) * 2 * ps * (1 - ps), n * (1 - ps)^2 + n*f*ps*(1-ps))
         Xt <- rbind(Xt, X)
       }
    }
    else {# sampling conditional on allele frequency
       if(is.null(p)) p <- 0.5
       nAll <- 2*n
       nA <- round(p*nAll,digits=0)
       nB <- nAll - nA
       Pop <- c(rep(1,nA),rep(0,nB))
       Xt <- NULL
       for (i in 1:nm) {
          sam <- matrix(sample(Pop,nAll),ncol=2,byrow=TRUE)
          status <- apply(sam,1,sum)
          nAA <- sum(status==2)
          nAB <- sum(status==1)
          nBB <- sum(status==0)
          if ((nAA+nAB+nBB)!=n) stop("HWData: error")
          Xt <- rbind(Xt,c(nAA,nAB,nBB))
       }
    }
    Xc <- Xt/n
    colnames(Xt) <- c("AA", "AB", "BB")
    colnames(Xc) <- c("AA", "AB", "BB")
    return(list(Xt = Xt, Xc = Xc))
}

