HWExact <-
function(X,alternative="two.sided",pvaluetype="dost",verbose=FALSE) {
n <- sum(X)
nA <- 2*X[1]+X[2]
nB <- 2*n-nA
MaxHet <- min(nA, nB)
if (MaxHet < 2) { # monomorphic marker, or only one outcome
   pval <- 1
   prob <- 1 
   pofthesample <- 1
   ind <- 1
#  ind <- match(X[2],seq(MaxHet%%2,MaxHet,2))
}
else {
   ind <- match(X[2],seq(MaxHet%%2,MaxHet,2))
   enAB <- nA*nB/(2*n-1) # Expected number of heterozygotes
   enAB <- round(enAB,digits=0)
   if ((enAB%%2) != (MaxHet%%2)) enAB <- enAB +1 # enAB and MaxHet should both be odd or even
   nAA <- (nA-enAB)/2
   nBB <- (nB-enAB)/2
   initialprob <- 1 # cheaper
   #initialprob <- HWCondProbAB(n, nA, enAB)$p
   AboveExp <- NULL
   BelowExp <- NULL
   if(enAB < MaxHet) AboveExp <- CompProbUp(nAA,nBB,enAB,initialprob,MaxHet)
   BelowExp <- CompProbDown(nAA,nBB,enAB,initialprob)
   prob <- c(rev(BelowExp),initialprob,AboveExp)
   prob <- prob/sum(prob)
 }
   Plow <- cumsum(prob)
   Phigh <- 1-c(0,Plow)
   Phigh <- Phigh[-length(Phigh)]
   Phwe <- pmin(1,2*Phigh,2*Plow) 
   pofthesample <- prob[ind] # cheaper
# pofthesample <- HWCondProbAB(n, nA, X[2])$p
# next lines for debugging
# RR <- round(cbind(seq((MaxHet%%2),MaxHet,2),prob,cumsum(prob),Phwe,Phigh,Plow),digits=6)
# colnames(RR) <- c("nAB","P(NAB=nAB|nA)","Cum","Phwe","Phigh","Plow")
# print(RR)
pval <- switch(alternative,
               greater = Phigh[ind],
               less = Plow[ind],
               two.sided = switch(pvaluetype,
                                  dost = Phwe[ind], 
                                  selome = sum(prob[prob <= pofthesample]),
                                  stop("invalid value for parameter pvaluetype")),
               stop("invalid value for parameter alternative"))
if(verbose) {
   D <- HWChisq(X)$D
   cat("Haldane's Exact test for Hardy-Weinberg equilibrium\n")
   cat("sample counts: nAA = ",X[1],"nAB = ",X[2],"nBB = ",X[3],"\n")
   stringtwosided <- paste("H0: HWE (D==0), H1: D <> 0 \nD = ",format(D,scientific=FALSE),"p = ",
                           format(pval,scientific=FALSE),"\n")
   stringgreater <- paste("H0: HWE (D==0), H1: D > 0 \nD = ",format(D,scientific=FALSE),"p = ",
                          format(pval,scientific=FALSE),"\n")
   stringless <- paste("H0: HWE (D==0), H1: D < 0 \nD = ",format(D,scientific=FALSE),"p = ",
                       format(pval,scientific=FALSE),"\n")
   toprint <- switch(alternative,
      two.sided = stringtwosided,
      greater =   stringgreater,
      less =      stringless)
   cat(toprint)
}
return(list(pval=pval,prob=prob,pofthesample=pofthesample))
}

