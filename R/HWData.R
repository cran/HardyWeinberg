HWData <- function (nm=100,n=rep(100,nm),f=rep(0,nm),p=runif(nm), pfixed = FALSE, 
    exactequilibrium = FALSE, pdist = "runif", ...) 
{
    Xt <- NULL
    if(length(p)==1) p <- rep(p,nm)
    if(length(n)==1) n <- rep(n,nm)
    if(length(f)==1) f <- rep(f,nm)
    ps <- p
    # allele frequency as a random variable
    if (!pfixed) {
        for (i in 1:nm) {
            if (is.null(p[i])) {
                if (pdist == "runif") 
                  ps[i] <- runif(1, ...)
                else {
                  if (pdist == "rbeta") 
                    ps[i] <- rbeta(1, ...)
                  else stop("unknown value for pdist")
                }
            }
            else ps <- p
            if (!exactequilibrium)
              # sample from multinomial distribution with HWE genotype frequencies.
                X <- t(rmultinom(1, size = n[i], prob = c(ps[i]^2 + 
                  f[i] * ps[i] * (1 - ps[i]), (1 - f[i]) * 2 * ps[i] * (1 - 
                  ps[i]), (1 - ps[i])^2 + f[i] * ps[i] * (1 - ps[i]))))
            else X <- c(n[i] * ps[i]^2 + n[i] * f[i] * ps[i] * (1 - ps[i]), n[i] * 
                (1 - f[i]) * 2 * ps[i] * (1 - ps[i]), n[i] * (1 - ps[i])^2 + 
                n[i] * f[i] * ps[i] * (1 - ps[i]))
            Xt <- rbind(Xt, X)
        }
    }
    else {
      # allele frequency is fixed
        if (is.null(p)) 
            p <- rep(0.5,nm)
        nAll <- 2 * n
        nA <- round(p * nAll, digits = 0)
        nB <- nAll - nA
        Xt <- NULL
        for (i in 1:nm) {
            Pop <- c(rep(1, nA[i]), rep(0, nB[i]))
            sam <- matrix(sample(Pop, nAll), ncol = 2, byrow = TRUE)
            status <- apply(sam, 1, sum)
            nAA <- sum(status == 2)
            nAB <- sum(status == 1)
            nBB <- sum(status == 0)
            if ((nAA + nAB + nBB) != n[i]) 
                stop("HWData: error")
            Xt <- rbind(Xt, c(nAA, nAB, nBB))
        }
    }
    Xc <- Xt/n
    colnames(Xt) <- c("AA", "AB", "BB")
    colnames(Xc) <- c("AA", "AB", "BB")
    return(list(Xt = Xt, Xc = Xc))
}
