HWEM <-
function(x,p=NULL,G=NULL,delta.init=c(0.5,0.5),itmax=50,
                 eps=1e-6,verbose=FALSE) {
  S <- length(delta.init)
  delta.init <- delta.init/sum(delta.init)
  
  if(is.null(p) & is.null(G)) {
    stop("Either allele frequencies (p) or genotype frequencies (G) need to be specified for each group.")
  }
  
  if(!is.null(p) & !is.vector(p)) {
    stop("Allele frequencies (p) should be given as a single scalar or a vector")
  }
  
  if(is.null(G)) {
    G <- matrix(NA,length(x),length(p))
    for(i in 1:length(p)) {
      G[,i] <- c(p[i]*p[i],2*p[i]*(1-p[i]),(1-p[i])*(1-p[i]))
    }
  }

  Dx <- diag(x/sum(x))
  
  change.ssq <- 100
  iter <- 1
  while(change.ssq > eps & iter < itmax) {
    Ddelta <- diag(delta.init)
    G1  <- G%*%Ddelta
    G1  <- G1/rowSums(G1)
    G1  <- Dx%*%G1
    delta.new <- colSums(G1)
    e <- delta.new-delta.init
    change.ssq <- sum(e*e)
    if(verbose) {
      cat (
          "iteration ",
          formatC (
            iter,
            digits = 0,
            width = 3,
            format = "d"
          ),
          " group1",
          formatC (
            delta.new[1],
            digits = 6,
            width = 10,
            format = "f"
          ),
          " group2",
          formatC (
            delta.new[2],
            digits = 6,
            width = 10,
            format = "f"
          ),
          " change.ssq",
          formatC (
            change.ssq,
            digits = 8,
            width = 10,
            format = "f"
          ),
          "\n"
        )
    }
    iter <- iter+1
    delta.init <- delta.new
  }
  names(delta.new) <- paste("Group",1:length(delta.new),sep="")
  return(delta.new)
}
