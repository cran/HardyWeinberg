MakeCounts <-function(X,alleles,pos1=1,pos2=3) {
  n <- nrow(X)
  p <- ncol(X)
  Y <- NULL
  for(j in 1:p) {
    snp <- X[,j]
    al1 <- substr(alleles[j],pos1,pos1)
    al2 <- substr(alleles[j],pos2,pos2)
    homAA <- paste(al1,al1,sep="")
    homBB <- paste(al2,al2,sep="")
    hetAB <- paste(al1,al2,sep="")
    nAA <- sum(snp==homAA,na.rm=TRUE)
    nBB <- sum(snp==homBB,na.rm=TRUE)
    nAB <- sum(snp==hetAB,na.rm=TRUE)
    nNA <- sum(is.na(snp))
    tot <- nAA+nAB+nBB+nNA
    if(tot!=n) {
      cat(j,"\n")
      stop("genotypes and missings do not sum n")
    }
    Y <- rbind(Y,c(nAA,nAB,nBB,nNA))
  }
  colnames(Y) <- c("AA","AB","BB","NA")
  return(Y)
}
