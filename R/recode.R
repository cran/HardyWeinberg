recode <- function(X,alleles,values=c(0,1,2),pos1=1,pos2=3) {
   n <- nrow(X)
   p <- ncol(X)
   Y <- NULL
   for(i in 1:p) {
      snp <- rep(NA,n)
      # extract first allele
      al1 <- substr(alleles[i],pos1,pos1)
      # extract second allele
      al2 <- substr(alleles[i],pos2,pos2)
      if(i%%100 ==0) cat("Converting marker ",i,"\n")
      hom1 <- paste(al1,al1,sep="")
      hom2 <- paste(al2,al2,sep="")
      het <- paste(al1,al2,sep="")
      snp[X[,i]==hom1] <- values[1]
      snp[X[,i]==het] <- values[2]
      snp[X[,i]==hom2] <- values[3]
      Y <- cbind(Y,snp)
    }
   colnames(Y) <- paste("SNP",1:p,sep="")
   return(Y)
}
