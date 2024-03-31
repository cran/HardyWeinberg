HWf <- function(X) 
{
  if(is.vector(X)) {
    if (length(X) != 3 | any(X < 0)) 
      stop("HWf: X is not a 3 by 1 non-negative count vector")
    if (any(!is.wholenumber(X))) {
      warning("Genotype counts are not integers, counts will be rounded.")
      X <- round(X, digits = 0)
    }
    X <- order.auto(X)
    if (maf(X) == 0) {
      warning("Monomorphic marker, f is NaN.")
      fhat <- NaN
    }
    else {
      nA <- 2 * X[1] + X[2]
      nB <- 2 * X[3] + X[2]
      fhat <- (4 * X[1] * X[3] - (X[2])^2)/(nA * nB)
      names(fhat) <- NULL
    }
  } # end if(is.vector)
  if(is.matrix(X)) { # its a matrix of genotype counts
    X <- X[,names(order.auto(X[1,]))]
    nA <- 2 * X[, 1] + X[, 2]
    nB <- 2 * X[, 3] + X[, 2]
    fhat <- (4 * X[, 1] * X[, 3] - X[, 2] * X[, 2])/(nA * nB)
  }
  return(fhat)
}


