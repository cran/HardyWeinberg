mac <- function (X) 
{
  # compute the minor allele count for a vector or for each row of matrix of genotype counts.
  if(is.vector(X)) {
    X <- order.auto(X)
    y <- 2*X[1] + X[2] # number of minor alleles
  } else if (is.matrix(X) | is.data.frame(X)) {
    X <- X[,names(order.auto(X[1,]))]
    n <- apply(X,1,sum)
    nA <- 2*X[,1] + X[,2]
    nB <- 2*X[,3] + X[,2]
    y <- pmin(nA, nB)
  } else {
    stop("unknown input format.")
  }
  return(y)
}