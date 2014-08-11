vaf <- function(X,hw=FALSE,minor=TRUE) {
  n <- sum(X)
  if (length(X) != 3 | any(X < 0)) 
     stop("vaf: X is not a 3 by 1 non-negative count vector")
  if (any(!is.wholenumber(X))) {
     warning("Genotype counts are not integers, counts will be rounded.")
     X <- round(X, digits = 0)
  }
  Xhom <- X[homozyg(X)]
  Xhet <- X[heterozyg(X)]
  if(minor) X <- c(min(Xhom), Xhet, max(Xhom)) else X <- c(max(Xhom), Xhet, min(Xhom))
  names(X) <- c("AA", "AB", "BB")
  pA <- af(X)
  if (is.vector(X)) {
     pAA <- X[1]/sum(X)
  }
  else if (is.matrix(X)) {
     pAA <- X[,1]/apply(X, 1, sum)
  }
  if(!hw) y <- (pA+pAA-2*(pA^2))/(2*n) else y <- pA*(1-pA)/(2*n)
  return(y)
}
