is.mono <- function(x) {
  if(is.vector(x)) {
    if(length(x) == 3) { # biallelic autosomal
      x <- order.auto(x)
      n <- sum(x)
      nA <- 2*x[1] + x[2]
      nB <- 2*x[3] + x[2]
      y <- (nA == 2*n) | (nB == 2*n)
    } else if(length(x) == 5) { # biallelic X-chromosomal
      x <- order.x(x)
      nt <- sum(x[1:2]) + 2*sum(x[3:5])
      nA <- x[1] + 2*x[3] + x[4]
      nB <- x[2] + 2*x[5] + x[4]
      y <- (nA == nt) | (nB == nt)
    } else {
      stop("x should have length 3 or 5.")
    }
  } else if(is.matrix(x)) {
    if(ncol(x) == 3) { # biallelic autosomal
      x <- x[,names(order.auto(x[1,]))]
      n <- apply(x,1,sum)
      nA <- 2*x[,1] + x[,2]
      nB <- 2*x[,3] + x[,2]
      y <- (nA == 2*n) | (nB == 2*n)
    } else if(ncol(x) == 5) { # biallelic X-chromosomal
      x <- x[,names(order.x(x[1,]))]
      nt <- x[,1] + x[,2] + 2*rowSums(x[,3:5])
      nA <- x[,1] + 2*x[,3] + x[,4]
      nB <- x[,2] + 2*x[,5] + x[,4]
      y <- (nA == nt) | (nB == nt)
    } else {
      stop("x should have 3 or 5 columns.")
    }
  } else if(is.data.frame(x)) {
    x <- as.matrix(x)
    if(ncol(x) == 3) { # biallelic autosomal
      x <- x[,names(order.auto(x[1,]))]
      n <- apply(x,1,sum)
      nA <- 2*x[,1] + x[,2]
      nB <- 2*x[,3] + x[,2]
      y <- (nA == 2*n) | (nB == 2*n)
    } else if(ncol(x) == 5) { # biallelic X-chromosomal
      x <- x[,names(order.x(x[1,]))]
      nt <- x[,1] + x[,2] + 2*rowSums(x[,3:5])
      nA <- x[,1] + 2*x[,3] + x[,4]
      nB <- x[,2] + 2*x[,5] + x[,4]
      y <- (nA == nt) | (nB == nt)
    } else {
      stop("x should have 3 or 5 columns.")
    }
  } else {
    y <- NULL
    stop("is.mono: unknown input format; supply a matrix or a vector.")
  }
  names(y) <- NULL
  return(y)
}
