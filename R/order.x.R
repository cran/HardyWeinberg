order.x <- function(X) {
  theals <- alleles(X)
  if(length(theals) != 2) {
    stop("The marker is not bi-allelic X-chromosomal.\n Supply counts as a named vector c(A,B,AA,AB,BB)")
  }
  hom1 <- paste(theals[1],theals[1],sep="")
  hom2 <- paste(theals[2],theals[2],sep="")
  type1 <- c(theals[1],hom1)
  type2 <- c(theals[2],hom2)
  t1 <- X[type1]
  t2 <- X[type2]
  nt1 <- t1[1] + 2*t1[2]
  nt2 <- t2[1] + 2*t2[2]
  if(nt1 <= nt2) {
    minoral <- names(nt1)
    majoral <- names(nt2)
  } else {
    minoral <- names(nt2)
    majoral <- names(nt1)
  }
  minorhom <- paste(minoral,minoral,sep="")
  majorhom <- paste(majoral,majoral,sep="")
  diploids <- names(X)[nchar(names(X))==2]
  Y <- X[diploids]
  ind.het <- heterozyg(Y)
  thehet <- names(Y)[ind.het]
  newseq <- c(minoral,majoral,minorhom,thehet,majorhom)
  Z <- X[newseq]
  return(Z)
}
