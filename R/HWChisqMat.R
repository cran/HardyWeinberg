HWChisqMat <-
function(X,...) {

pvalvec <- NULL
chisqvec <- NULL
Dvec <- NULL

for (i in 1:nrow(X)) {
   out <- HWChisq(X[i,],...)
   pvalvec <- c(pvalvec,out$pval)
   chisqvec <- c(chisqvec,out$chisq)
   Dvec <- c(Dvec,out$D)   
}
return(list(pvalvec=pvalvec,chisqvec=chisqvec,Dvec=Dvec))
}

