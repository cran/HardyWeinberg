`HWData` <-
function(n=100,nm=100,pfixed=NULL,exactequilibrium=FALSE,pdist="runif",...) {
#
# Generates data in HWE
#
Xt <- NULL
for (i in 1:nm) {
   if(is.null(pfixed)) {
     if(pdist=="runif") p <- runif(1,...)
     else {
       if(pdist=="rbeta")
         p <- rbeta(1,...)
       else
         stop("unknown value for pdist")
     }
   }
   else p <- pfixed
   if(!exactequilibrium)
     X <- t(rmultinom(1, size = n, prob=c(p^2,2*p*(1-p),(1-p)^2))) else
     X <- c(n*p^2,n*2*p*(1-p),n*(1-p)^2)  
   Xt <- rbind(Xt,X)
}
Xc <- Xt/n
colnames(Xt) <- c("AA","AB","BB")
colnames(Xc) <- c("AA","AB","BB")
return(list(Xt=Xt,Xc=Xc))
}

