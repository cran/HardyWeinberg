mac <- function (X) 
{
    # compute the minor allele count for a vector or for each row of matrix of genotype counts.
    Xhom <- X[homozyg(X)]
    Xhet <- X[heterozyg(X)]
    X <- c(min(Xhom), Xhet, max(Xhom))
    if (is.vector(X)) {
        n <- sum(X)
        nA <- 2*X[1] + X[2] # number of A alleles
        nB <- 2*n - nA
        y <- min(nA,nB)
    }
    else if (is.matrix(X)) {
        n <- apply(X,1,sum)
        nA <- 2*X[,1] + X[,2]
        nB <- 2*n - nA
        y <- pmin(nA,nB)
    }
    return(y)
}
