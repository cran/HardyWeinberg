`HWExact` <-
function(X,singleterms=FALSE,alternative="two.sided",sorting="probability",verbose=FALSE) {
  nAAs <- X[1]
  nABs <- X[2]
  nBBs <- X[3]
  n <- sum(X)
  nAs <- 2*nAAs+nABs
  nBs <- 2*nBBs+nABs
  MaxHet <- min(nAs,nBs)
  D <- HWChisq(X)$D
  pofthesample <- HWCondProbAB(n,nAs,nABs)$p

  if(MaxHet%%2==0) nABvec <- seq(0,MaxHet,by=2)
  if(MaxHet%%2==1) nABvec <- seq(1,MaxHet,by=2)

  singterm <- NULL
  cumterm <- NULL
  nAAvec <- NULL
  nBBvec <- NULL
  Chivec <- NULL
  pvalvec <- NULL
  Chiccvec <- NULL
  pvalccvec <- NULL
  Dvec <- NULL
     
  for (i in 1:length(nABvec)) {
     out <- HWCondProbAB(n,nAs,nABvec[i])
     p <- out$p
     singterm <- c(singterm,p)
     nAB <- nABvec[i]
     nAA <- 0.5*(nAs-nABvec[i])
     nBB <- 0.5*(nBs-nABvec[i])
     nAAvec <- c(nAAvec,nAA)
     nBBvec <- c(nBBvec,nBB)
     samp <- c(nAA,nAB,nBB)
     out <- HWChisq(samp)
     Dvec <- c(Dvec,out$D)
     Chivec <- c(Chivec,out$chisq)
     pvalvec <- c(pvalvec,out$pval)
  }

  if(alternative=="two.sided") {
     pval <- sum(singterm[singterm<=pofthesample])
  }
  if(alternative!="two.sided") {
     stop("unkown option for alternative")
  }
  
  if(verbose) {
     cat("Exact test for Hardy-Weinberg equilibrium\n")
     cat("sample counts: nAA = ",nAAs,"nAB = ",nABs,"nBB = ",nBBs,"\n")
  
     if(alternative=="two.sided") {
         cat("H0: HWE (D=0), H1: D <> 0 \n")
         cat("D = ",D,"p = ",format(pval,scientific=FALSE),"\n")
     }
  }

  if(sorting=="probability") ind <- order(singterm)
  else ind <- order(nABvec)

  cumterm <- cumsum(singterm)

  if(singleterms) {

     Result <- data.frame(nAAvec,nABvec,nBBvec,round(cbind(singterm,cumterm,Chivec,pvalvec,Dvec),digits=8),
                          row.names=1:length(nABvec))
     colnames(Result) <- c("AA","AB","BB","Single term","Prob","X2","pval","D")

     if(verbose) {
        Result <- Result[ind,]
        Result$Prob <- cumsum(Result$"Single term")
        cat("\nProbabilities and statistics for all possible samples:\n")
        print(Result)
      }
  }

  return(list(D=D,pval=pval,pofthesample=pofthesample))  
}

