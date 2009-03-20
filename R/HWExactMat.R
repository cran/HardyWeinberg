HWExactMat <-
function(X,...) {

pvalvec <- NULL
nmarkers <- nrow(X)

for (i in 1:nrow(X)) {
   out <- HWExact(X[i,],...)
   pvalvec <- c(pvalvec,out$pval)
}

return(list(pvalvec=pvalvec))
}

