MakeColourVector <-
function(x) {
  #
  # create a vector of colours according to a third quantitative 
  # variable; this variable should be in the (0,1) interval
  #
  cf <- colorRamp(c("green","yellow","red"))

  M <- cf(x)/255
  n <- nrow(M)
  colvec <- numeric(n)
  for(i in 1:n) {
    colvec[i] <- rgb(M[i,1],M[i,2],M[i,3])
  }
  return(colvec)
}
