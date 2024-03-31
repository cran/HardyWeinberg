he <-
function(x,bias.correct=TRUE) {
  # calculate expected heterozygosity
  if(is.vector(x)) {
    n <- sum(x)
    p <- af(x)
    q <- 1 - p
    eh <- 2*p*q
    if(bias.correct) {
      eh <- 2*n*eh/(2*n-1)
    }
    names(eh) <- "He"  
  } else if(is.matrix(x)) {
    n <- rowSums(x)
    p <- af(x)
    q <- 1 - p
    eh <- 2*p*q
    if(bias.correct) {
      eh <- 2*n*eh/(2*n-1)
    }
  }
  return(eh)
}
