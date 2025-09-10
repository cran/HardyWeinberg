HWBalance <- function(x,verbose=TRUE,alpha=0.05,alternative="two.sided") {
  level <- 100*(1-alpha)
  values <- c("two.sided","greater","less","excess","lack")
  if(!(alternative %in% values)) {
    stop("improper value for alternative,\n")
  }
  if(is.vector(x) & length(x)==3) { # biallelic case
    if(any(x==0)) {
      cat("This test does not allow for zero genotype counts.\n")
      cat("Replace zeros or use a different test.\n")
    }
    n <- sum(x)
    g <- x/n
    K <- 2
    Z1 <- sqrt((K*(K-1))/(K+1))*(log(g[2]) - 0.5*(log(g[1]) + log(g[3])))
    VZ1 <- (K-1)*(1/g[1] + 1/g[3]) + (4/(K-1))*(1/g[2])
    VZ1 <- VZ1/(n*K*(K+1))
    refvalue <- sqrt(K*(K-1)/(K+1))*log(2)
    Z <- (Z1 - refvalue)/sqrt(VZ1)
    if(alternative=="two.sided") {
      pval <- 2*pnorm(abs(Z),lower.tail = FALSE)
      ll <- Z1 - qnorm(1-alpha/2)*sqrt(VZ1)
      ul <- Z1 + qnorm(1-alpha/2)*sqrt(VZ1)
    } else if (alternative=="greater") {
      pval <- pnorm(Z,lower.tail = FALSE)
      ll <- Z1 - qnorm(1-alpha)*sqrt(VZ1)
      ul <- Inf
    } else if (alternative=="less") {
      pval <- pnorm(Z,lower.tail = TRUE)
      ll <- -Inf
      ul <- Z1 + qnorm(1-alpha)*sqrt(VZ1)
    }
    chisq <- Z*Z
  } else if (is.matrix(x)) {
      n <- sum(x)
      K <- ncol(x)
      refvalue <- sqrt(K*(K-1)/(K+1))*log(2)
      if(any(x[lower.tri(x)]==0)) {
        cat("This test does not allow for zero genotype counts.\n")
        cat("Replace zeros or use a different test.\n")
      }
      GF <- x/sum(x)
      hom.gf <- diag(GF)
      het.gf <- GF[lower.tri(GF)]
      sum.recip.hom <- sum(1/hom.gf)
      sum.recip.het <- sum(1/het.gf)
      numvar <- (K-1)*sum.recip.hom + (4*sum.recip.het)/(K-1)
      VZ1 <- numvar/(n*K*(K+1))
      log.geom.hom <- mean(log(hom.gf))
      log.geom.het <- mean(log(het.gf))
      Z1 <- sqrt(K*(K-1)/(K+1))*(log.geom.het - log.geom.hom)
      Z <- (Z1 - refvalue)/sqrt(VZ1)
      if(alternative=="two.sided") {
         pval <- 2*pnorm(abs(Z),lower.tail = FALSE)
         ll <- Z1 - qnorm(1-alpha/2)*sqrt(VZ1)
         ul <- Z1 + qnorm(1-alpha/2)*sqrt(VZ1)
      } else if (alternative=="greater") {
         pval <- pnorm(Z,lower.tail = FALSE)
         ll <- Z1 - qnorm(1-alpha)*sqrt(VZ1)
         ul <- Inf
      } else if (alternative=="less") {
         pval <- pnorm(Z,lower.tail = TRUE)
         ll <- -Inf
         ul <- Z1 + qnorm(1-alpha)*sqrt(VZ1)
      }
      chisq <- Z*Z
      
  } else {
      cat("Inappropriate data format. Supply x as 3 x 1 vector or a K x K lower triangular matrix of genotype counts.\n")
  }
  stringtwosided <- "H0: HWE (z_1==reference value), H1: z_1 <> reference value \n"
  stringgreater  <- "H0: HWE (z_1==reference value), H1: z_1 > reference value \n" 
  stringless     <- "H0: HWE (z_1==reference value), H1: z_1 < reference value \n" 
  if(verbose) {
      cat("Asymptotic compositional test for HWE based on the heterozygote-homozygote balance.\n")
      cat(K,"alleles detected.\n")
      cat("Equilibrium reference value:",round(refvalue,4),"\n")
      cat("Sample heterozygote-homozygote balance \\hat{z}_1 = ",Z1,"\n")
      toprint <- switch(alternative, two.sided = stringtwosided, 
                greater = stringgreater, less = stringless)
      cat(toprint)
      cat("Z = ",Z,"CHISQ = ",chisq,"p-value = ",pval,"\n")
      cat(level,"% CI(z_1) = (",round(ll,4),",",round(ul,4),")\n")
    }
  out <- list(Z1=Z1,VZ1=VZ1,Z=Z,chisq=chisq,pval=pval)
}
