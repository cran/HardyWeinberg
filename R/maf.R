maf <- function(x,option=1,verbose=FALSE) {
  if(is.vector(x)) {
    layout <- 1 # a single marker
  } else if (is.matrix(x)) {
    if(nrow(x) == ncol(x)) {
      layout <- 2 # a single marker in matrix format
    } else {
      layout <- 3 # multiple markers in rows
    }
  } else {
    stop("Unkown format, x must be a matrix or a vector")
  }
  if(verbose) {
    if(layout==1) cat("input is a vector of genotype counts for a single marker.\n")
    else if(layout==2) cat("input is a matrix of genotype counts for a single marker.\n")
    else if(layout==3) cat("input is a matrix with 3 columns of genotype counts for a set of markers, one per row.\n")
  }
  if(layout==1) { # a single marker in a vector
    if(is.null(names(x))) {
      warning("Genotype counts are not labelled, default sequence (AA,AB,BB) assumed.")
      x <- genlabels(x)
    }
    n   <- sum(x,na.rm=TRUE)
    na  <- n.alleles(x)
    als <- alleles(x)
    if(verbose) {
      cat(na,"alleles detected\n")
      print(als)
    }
    genot <- expand.grid(alleles(x),alleles(x))
    genot <- paste(genot$Var2,genot$Var1,sep="")
    GT <- matrix(genot,ncol=na,byrow = TRUE) # names of genotypes
    rownames(GT) <- als
    colnames(GT) <- rownames(GT)
    GTc <- matrix(0,nrow=nrow(GT),ncol=ncol(GT)) # counts of genotypes
    rownames(GTc) <- rownames(GT)
    colnames(GTc) <- colnames(GT)
    for(i in 1:na) {
      for(j in 1:na) {
        ind.gt <- names(x)==GT[i,j]
        GTc[i,j] <- sum(x[ind.gt],na.rm=TRUE)
      }
    }
    if(verbose) {
      cat("genotype frequencies:\n")
      print(GTc)
    }
    acounts <- numeric(na)
    names(acounts) <- als
    for(i in 1:na) {
      acounts[i] <- sum(GTc[i,] + GTc[,i])  
    }
    acounts <- sort(acounts)
    nt <- as.integer(sum(acounts))
    if(nt != as.integer(2*n)) {
      stop("number of alleles does not double the sample size")
    }
  } else if(layout==2) { # it's square matrix of genotype counts.
    n  <- sum(x)
    na <- nrow(x)
    acounts <- numeric(na)
    acounts <- rowSums(x) + colSums(x)
    acounts <- sort(acounts)
    nt <- as.integer(sum(acounts))
    twon <- as.integer(2*n)
    if(twon!=nt) stop("number of alleles is not twice the sample size.")
  } else if(layout==3) { # it's a matrix with multiple biallelic markers; 
    if(is.null(colnames(x))) {
      warning("Genotype counts are not labelled, default labels (AA,AB,BB) assumed.")
      x <- genlabels(x)
    } else { # if colnames given, order the columns
      x <- x[,names(order.auto(x[1,]))]
    }
    n <- rowSums(x) # sample sizes may differ
    acounts <- mac(x)
    names(acounts) <- rownames(x)
  } else {
    stop("unknown layout.")
  }
  if(option == 1) {
    if(layout != 3) {
      y <- acounts[which.min(acounts)]/(2*n)  
    } else {
      y <- acounts/(2*n)
    }
  } else if (option == 2) {
    if(layout != 3) {
      y <- acounts/(2*n)
    } else {
      y <- acounts/(2*n)
      y <- cbind(y,1-y)
      colnames(y) <- c("Minor","Major")
    }
  } else if (option == 3) {
    if(layout !=3) {
      y <- acounts
    } else {
      y <- cbind(acounts,2*n-acounts)
      colnames(y) <- c("Minor","Major")
    }
  } else {
    stop("unknown parameter option")
  }
  return(y)
}

