order.auto <- function(X) {
  if(is.null(names(X))) {
    warning("No genotype labels supplied; default order (AA,AB,BB) is assumed.")
    names(X) <- c("AA","AB","BB") 
  }
  thehet   <- names(X)[heterozyg(X)]  
  homs     <- sort(X[homozyg(X)])
  minorhom <- names(homs)[1]
  majorhom <- names(homs)[2]
  newseq   <- c(minorhom, thehet, majorhom)
  return(X[newseq])
}
