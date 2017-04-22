HWExactStats <- function(X,x.linked=FALSE,plinkcode=TRUE,midp=FALSE,...) {
  X <- as.matrix(X)
  nmarkers <- nrow(X)
  pvalvec <- numeric(nmarkers)
  if(plinkcode & !x.linked) {
    for (i in 1:nmarkers) {
      pvalvec[i] <- SNPHWE2(X[i,2],X[i,1],X[i,3],midp)
    }
  } else {
    for (i in 1:nmarkers) {
      pvalvec[i] <- HWExact(X[i,],x.linked=x.linked,...)$pval
    }  
  }
  return(pvalvec)
}
