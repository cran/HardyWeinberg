af <- function(x) {
  if(is.vector(x)) {
    lx <- length(x)
    if(lx != 3) {
      stop("x does not contain frequencies or counts of a biallelic autosomal marker")
    } 
    if(is.null(names(x))) x <- genlabels(x)
    if(!heterozyg(x)[2]) cat("Warning: the second count is not a heterozygote count.\n") 
    p <- (x[1]+0.5*x[2])/sum(x)  
    names(p) <- substr(names(x)[1],1,1)
  } else if(is.matrix(x)) {
    nc <- ncol(x)
    if(nc != 3) {
      stop("x does not contain frequencies or counts of biallelic autosomal markers")
    }
    if(is.null(colnames(x))) x <- genlabels(x)
    if(!heterozyg(x[1,])[2]) cat("Warning: the second column is not a heterozygote count.\n") 
    p <- (x[,1]+0.5*x[,2])/apply(x,1,sum)
  }
  return(p)
}
