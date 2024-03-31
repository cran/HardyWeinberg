HWLindley.cri <- function(x,verbose=TRUE,limits=c(0.025,0.975)) {
  #
  # calculate a credible interval, by default 95%, for Lindley's posterior.
  #
  theint <- function(x,alpha=0,limits=c(0.025,0.975)) {
    z <- integrate(HWLindley,-Inf,alpha,x)
    value <- z$value
    pro <- value
    return(pro)
  }
  Lindley.pct <- function(alpha,x,value) {
    y <- theint(x,alpha) - value
  }
  x <- order.auto(x)
  ll <- uniroot(Lindley.pct,c(-2,2),x,limits[1])
  ll <- ll$root
  ul <- uniroot(Lindley.pct,c(-2,2),x,limits[2])
  ul <- ul$root
  outputstring <- paste("(",toString(round(ll,4)),",",
                            toString(round(ul,4)),")",sep="")
  percent <- round(100*(limits[2] - limits[1]),digits=0)
  if(verbose) {
    cat(paste("Lindley's ",toString(percent),"% credible interval\n",sep=""))
    cat("HWE parameter alpha: ",outputstring,"\n")
  }
  out <- c(ll,ul)
}
