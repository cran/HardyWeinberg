HWNetwork <- function(a1,a2,ma=NULL,fe=NULL,gender=NULL,verbose=TRUE) {
  if(is.null(ma) & is.null(fe)) { # format of original genotype data
    
    alleles <- unique(sort(c(a1,a2)))
    k <- length(alleles) # the number of alleles
    
    allelenames <- paste("A",as.character(alleles),sep="")
    acounts <- numeric(k)
    for(i in 1:length(alleles)) {
      acounts[i] <- sum(a1==alleles[i],na.rm=TRUE) + 
        sum(a2==alleles[i],na.rm=TRUE)
    }
    names(acounts) <- allelenames
    
    ii <- order(acounts,decreasing=TRUE)
    acounts <- acounts[ii]
  
    alleles <- alleles[ii]
    allelenames <- allelenames[ii]
    
    acounts <- as.integer(acounts)
    names(acounts) <- allelenames
    
    if(is.null(gender)) {
      # autosomal test
    } else {
      # X-chromosomal test.
      n <- length(gender)
      if(any(gender!=1 & gender!=2)) stop("gender not properly coded (1 for males, 2 for females)")
      nm <- sum(gender==1,na.rm=TRUE)
      nf <- sum(gender==2,na.rm=TRUE)
     
      MaleFemaleCounts <- c(nm,nf)
      Pvalues <- c(0,0,0,0)
      
      acounts.males <- numeric(k)
      for(i in 1:length(alleles)) {
        acounts.males[i] <- sum(a1[gender==1]==alleles[i],na.rm=TRUE) + 
          sum(a2[gender==1]==alleles[i],na.rm=TRUE)
      }
      names(acounts.males) <- allelenames
      
      acounts.females <- numeric(k)
      for(i in 1:length(alleles)) {
        acounts.females[i] <- sum(a1[gender==2]==alleles[i],na.rm=TRUE) + 
          sum(a2[gender==2]==alleles[i],na.rm=TRUE)
      }
      names(acounts.females) <- allelenames
      
      #
      # Female genotype counts
      #
      
      f.a1 <- paste("A",as.character(a1[gender==2]),sep="")
      f.a2 <- paste("A",as.character(a2[gender==2]),sep="")
      fem <- paste(f.a1,f.a2,sep="/")
      
      fa1 <- factor(f.a1,levels=allelenames)
      fa2 <- factor(f.a2,levels=allelenames)
      
      M <- table(fa1,fa2)
      Mn <- as.matrix(unclass(M))
      Fec <- fold(Mn)
      
      prob.of.sample <- density.ma.gender(acounts.males,Fec)
      
      
      ostats <- c(0, prob.of.sample, 0, 0);  # input 4: prob (important!)
      observed=as.double(ostats)
      
      
      pval <- xChromosomal(acounts, MaleFemaleCounts, 
                           k, observed, 
                           Pvalues, 0, 0, 0, 0, 0, 0)
    }
      
  } else { # summarised data as tables of genotype counts
    
    nm <- sum(ma,na.rm=TRUE)
    nf <- sum(fe,na.rm=TRUE)
    
    MaleFemaleCounts <- c(nm,nf)
    
    k <- length(ma)
    
    Pvalues <- c(0,0,0,0)
    
    ft <- rowSums(fe) + colSums(fe) # the way to calculate female allele counts from
    # the triangular table; homozygotes are counted twice.
    
    acounts <- ma + ft
    
    acounts.males <- ma
    acounts.females <- ft
    
    ii <- order(acounts,decreasing = TRUE)
    
    acounts <- acounts[ii]
    acounts.males <- acounts.males[ii]
    acounts.females <- acounts.females[ii]
    
    prob.of.sample <- density.ma.gender(ma,fe)
    
    ostats <- c(0, prob.of.sample, 0, 0);  # input 4: prob (important!)
    observed=as.double(ostats)
    
    acounts <- as.integer(acounts)

    pval <- xChromosomal(acounts, MaleFemaleCounts, 
                         k, observed, 
                         Pvalues, 0, 0, 0, 0, 0, 0)
    
  }
  if(verbose) {
   cat("Network algorithm for HWE Exact test with multiple alleles\n")
   cat(k,"alleles detected.\n")
   cat(nm,"males and ",nf,"females\n")
   cat("Allele counts:\n")
   res.counts <- rbind(acounts.males,acounts.females,acounts)
   rownames(res.counts) <- c("Males","Females","All")
   print(res.counts)
   cat("Probability of the sample:",prob.of.sample,"\n")
   cat("p-value:",pval,"\n")
  }
  out <- list(pval=pval)
}
