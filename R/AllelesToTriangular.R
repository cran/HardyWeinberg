AllelesToTriangular <- function(A1,A2=NULL,given=NULL) {
  #
  # eliminate missings
  #
  if(is.null(A2)) {
    ii <- is.na(A1)
    A1 <- A1[!ii]
  } else {
    ii <- (is.na(A1) | is.na(A2))
    A1 <- A1[!ii]
    A2 <- A2[!ii]
  }
  if(is.null(A2)) { # one vector, so split
    nn <- length(A1)
    i1 <- seq(1,nn,by=2)
    i2 <- i1+1
    A1n <- A1[i1]
    A2n <- A1[i2]
  } else {
    A1n <- A1
    A2n <- A2
  }
  #
  # number of alleles
  #
  als <- unique(c(A1n,A2n))
  if(is.null(given)) {
    nals <- sort(paste("A",as.character(als),sep=""))  
  } else {
    nals <- given
  }
  A1n <- paste("A",A1n,sep="")
  A2n <- paste("A",A2n,sep="")
  k <- length(nals)
  nG <- length(A1n) # number of genotypes
  gen <- rep(NA,nG)
  for(i in 1:length(A1n)) {
    gen[i] <- paste(sort(c(A1n[i],A2n[i])),collapse ="/")
  }
  #
  # make genotypes
  #
  if(length(A1n)!=length(A2n)) {
    stop("Allele vectors have unequal length.")
  }
  M <- matrix(0,nrow=k,ncol=k)
  colnames(M) <- nals
  rownames(M) <- nals
  GTM <- M
  for(i in 1:k) {
    for(j in 1:(i)) {
      GTM[i,j] <- paste(sort(c(nals[i],nals[j])),collapse="/")
    }
  }
  for(i in 1:k) {
    for(j in 1:(i)) {
      M[i,j] <- sum(gen==GTM[i,j])
    }
  }
  return(M)
}
