`HWChisq` <-
function(X,cc = 0,verbose = FALSE) {

mono <- FALSE
cc <- rep(cc,3)
n <- sum(X)
p <- (2*X[1]+X[2])/(2*n)
q <- 1-p
if((p==0) | (q==0)) {
  if(verbose) cat("warning: monomorphic marker\n")
  mono <- TRUE
}
obs <- X
exp <- c(n*p^2,n*2*p*q,n*q^2)
D <- 0.5*(obs[2]-exp[2])
chi <- (abs(obs-exp) - cc)^2
if(!mono) {
   chi2 <- chi/exp
   chisq <- sum(chi2)
   pval <- 1-pchisq(chisq,1)
} else {
   chisq <- NA
   pval <- 1
}
if (verbose)
   cat("Chi2 = ",chisq,"p-value = ",pval,"D = ",D,"\n")
return(list(chisq=chisq,pval=pval,D=D,p=p))
}

