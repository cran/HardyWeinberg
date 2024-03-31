is.mono <- function(x) {
  if(is.vector(x)) {
    x <- order.auto(x)
    n <- sum(x)
    nA <- 2*x[1] + x[2]
    nB <- 2*x[3] + x[2]
    y <- (nA == 2*n) | (nB == 2*n)
  } else if(is.matrix(x)) {
    n <- apply(x,1,sum)
    x <- x[,names(order.auto(x[1,]))]
    nA <- 2*x[,1] + x[,2]
    nB <- 2*x[,3] + x[,2]
    y <- (nA == 2*n) | (nB == 2*n)
  } else if(is.data.frame(x)) {
    x <- as.matrix(x)
    n <- apply(x,1,sum)
    x <- x[,names(order.auto(x[1,]))]
    nA <- 2*x[,1] + x[,2]
    nB <- 2*x[,3] + x[,2]
    y <- (nA == 2*n) | (nB == 2*n)
  } else {
    y <- NULL
    stop("is.mono: unknown input format; supply a matrix or a vector.")
  }
  names(y) <- NULL
  return(y)
}
