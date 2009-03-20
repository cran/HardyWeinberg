`CritSam` <-
function(n=5,Dpos=TRUE,alphalimit=0.05) {
   X <- GenerateSamples(n)
   ncomp <- nrow(X)
   Res <- NULL
   Ds <- NULL
   pval <- NULL
   fA <- NULL
   for (i in 1:nrow(X)) {
#      m <- matrix(c(X[i,1],X[i,2]/2,X[i,2]/2,X[i,3]),ncol=2)
      fA <- c(fA,(2*X[i,1]+X[i,2])/(2*n))
      Ds <- c(Ds,HWChisq(X[i,])$D)
#      pval <- c(pval,fisher.test(m,alternative="two.sided")$p.value)
      pval <- c(pval,HWExact(X[i,],alternative="two.sided")$pval)
   }   

   Y <- data.frame(X[,1],X[,2],X[,3],fA,Ds,pval)
   colnames(Y) <- c("AA","AB","BB","fA","Ds","pval")
   if(Dpos) Y <- Y[Y$Ds>0,] else Y <- Y[Y$Ds<0,]
   Y <- Y[Y$pval<alphalimit,]
   fre <- unique(fA)
   for (i in 1:length(fre)) {
      Ys <- Y[Y$fA == fre[i],]
      if (nrow(Ys) > 0) {
         indi <- which.max(Ys$pval)
         Ys <- Ys[indi,]
         Res <- rbind(Res,c(Ys$AA,Ys$AB,Ys$BB))
      }
   }
   Xn <- Res/n
   return(list(Xn=Xn,Ds=Ds,fA=fA))
}

