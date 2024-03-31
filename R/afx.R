afx <- function(x) {
  if(is.vector(x)) {
    lx <- length(x)
    if(lx != 5) {
      stop("x does not contain frequencies or counts of a biallelic X-chromosomal marker")
    } 
    if(is.null(names(x))) x <- genlabels(x)
    if(!heterozyg(x[3:5])[2]) cat("Warning: the fourth count is not a heterozygote count.\n")
    nt <- x[1] + x[2] + 2*sum(x[3:5])
    p  <- (x[1] + 2*x[3] + x[4])/nt  
    names(p) <- substr(names(x)[1],1,1)
  } else if(is.matrix(x)) {
    nc <- ncol(x)
    if(nc != 5) {
      stop("x does not contain frequencies or counts of biallelic X-chromosomal markers")
    }
    if(is.null(colnames(x))) x <- genlabels(x)
    if(!heterozyg(x[1,3:5])[2]) cat("Warning: the fourth column is not a heterozygote count.\n")
    nt <- x[,1] + x[,2] + 2*(x[,3] + x[,4] + x[,5])
    p <- (x[,1] + 2*x[,3] + x[,4])/nt
  }
  return(p)
}
