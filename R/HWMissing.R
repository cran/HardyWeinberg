HWMissing <- function(X,imputecolumn=1,m=50,verbose=FALSE,fisher=FALSE,alpha=0.05,varest=1,...) {
  #
  # function for inference for HWP using Multiple imputation.
  #
     if (alpha <= 0 | alpha >= 1) 
       stop("HWMissing: alpha should be in the range (0,1)")
     objecttype <- class(X)
     choice <- switch(objecttype, data.frame = 1, factor = 2, stop("X should be a dataframe or a factor"))
     if(choice == 1) { 
        y  <- as.numeric(X[,imputecolumn]) # internally everything is recoded as AA, AB, BB.
        yf <- factor(y,levels=c(1,2,3),labels=c("AA","AB","BB"))
        X[,imputecolumn] <- yf
        imp <- mice(X,m=m,predictorMatrix=quickpred(X,minpuc=0.10),printFlag=FALSE,...)
     }
     if(choice == 2) { # impute assuming MCAR
        y <- as.numeric(X)
        yf <- factor(y,levels=c(1,2,3),labels=c("AA","AB","BB"))
        X <- data.frame(yf,rep(1,length(yf)))
        imp <- mice(X,m=m,predictorMatrix=quickpred(X,minpuc=0.10),printFlag=FALSE,method=c("sample",""),...)
     }
     Xmat <- NULL
     for(i in 1:m) {
        Ximp <- complete(imp,i)
        nAA <- sum(Ximp[,imputecolumn]=="AA")
        nAB <- sum(Ximp[,imputecolumn]=="AB")
        nBB <- sum(Ximp[,imputecolumn]=="BB")
        Ximp <- c(AA=nAA,AB=nAB,BB=nBB)
        Xmat <- rbind(Xmat,Ximp)
     }
     # combine estimates of f, using Rubin's pooling rules.
     cout <- combineC(Xmat,fisher=fisher,alpha=alpha,varest=varest)
     Res <- c(cout$fhatimp,cout$llf,cout$ulf,cout$pvalimp,cout$r,cout$gamma)
     names(Res) <- c("f","llci","ulci","p-value","r","gamma")
     if(verbose) {
        cat("Test for Hardy-Weinberg equilibrium in the presence of missing values\n")
        cat("Inbreeding coefficient f = ",round(cout$fhatimp,digits=4),"\n")
        cat(round(100*(1-alpha),digits=0),"% Confidence interval (",round(cout$llf,digits=4),",",round(cout$ulf,digits=4),")\n")
        cat("p-value = ",round(cout$pvalimp,digits=4),"\n")
        cat("Relative increase in variance of f due to missings: r = ",round(cout$r,digits=4),"\n")
        cat("Fraction of missing information about f: lambda = ",round(cout$gamma,digits=4),"\n")
     }
return(Res)
}

