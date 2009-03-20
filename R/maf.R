`maf` <-
function(x) {
  if(is.vector(x)) {
   p <- (x[1]+0.5*x[2])/sum(x)
   y <- min(p,1-p)
  } else
  if(is.matrix(x)) {
   p <- (x[,1]+0.5*x[,2])/apply(x,1,sum)
   y <- pmin(p,1-p)
  }
  return(y)
}

