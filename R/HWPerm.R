HWPerm <-
function(x, nperm=17000, verbose=TRUE, FUN=Chisquare, ...) {
  n <- sum(x)
  nA <- 2*x[1]+x[2]
  nB <- 2*n-nA
  stat.obs <- FUN(x)
  nlarger <- 0
  pseudodist <- NULL
  for(i in 1:nperm) {
     xx <- sample(c(rep("A",nA),rep("B",nB)))
     i1 <- seq(1,2*n,2)
     i2 <- seq(2,2*n,2)
     A1 <- xx[i1]
     A2 <- xx[i2]
     Geno <- paste(A1,A2,sep="")
     Geno[Geno=="BA"] <- "AB"
     nAA <- sum(Geno=="AA")
     nAB <- sum(Geno=="AB")
     nBB <- sum(Geno=="BB")
     y <- c(nAA,nAB,nBB)
     names(y) <- c("AA","AB","BB")
     stat.pseudo <- FUN(y)
     if(stat.pseudo >= stat.obs) nlarger <- nlarger + 1
     pseudodist <- c(pseudodist,stat.pseudo)
   }
   pval <- nlarger/nperm
   if(verbose) {
     cat("Permutation test for Hardy-Weinberg equilibrium\n")
     cat(nperm,"permutations. p-value:", pval,"\n")
   }
   return(pval)
}
