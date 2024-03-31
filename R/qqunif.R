qqunif <- function (x, logplot = FALSE, lbs = 1:length(x), texton = FALSE, 
          xylim = NULL, main = "Q-Q plot for a uniform distribution",
          plotline = 0, xlab = "Expected p-value", 
          ylab = "Observed p-value", colvec=rep("black",length(x)),
	  colline = "black", ...) 
{
  nmis <- sum(is.na(x))
  if(nmis>0) cat("qqunif:",nmis,"missing values excluded.\n")
  x <- x[!is.na(x)]
  xs <- sort(x, index.return = TRUE)
  colvec <- colvec[xs$ix]
  pvals <- xs$x
  n <- length(x)
  epvals <- (1:n - 0.5)/n
  if (!logplot) {
    plot(epvals, pvals, xlab = "Expected p-value", ylab = "Observed p-value", 
         xlim = c(0, 1), ylim = c(0, 1), main = main, col = colvec, ...)
    if (texton) 
      text(epvals, pvals, lbs)
    qp <- quantile(pvals, c(0.25, 0.75))
    qe <- quantile(epvals, c(0.25, 0.75))
    slope <- diff(qp)/diff(qe)
    int <- qp[1] - slope * qe[1]
    if(plotline==1) {
      abline(0, 1, lwd = 1, col = colline)  
    } else if(plotline==2) {
      abline(int, slope, lwd = 1, col = colline)
    } else {
      abline(int, slope, lwd = 1, col = colline)
    }
  }
  else {
    lpvals <- -log10(pvals)
    lepvals <- -log10(epvals)
    if (is.null(xylim)) 
      xylim <- max(c(lpvals, lepvals))
    plot(lepvals, lpvals, xlab = xlab, ylab = ylab,
         ylim = c(0, xylim), main = main, col = colvec, ...)
    if (texton) 
      text(lepvals, lpvals, lbs)
    qp <- quantile(lpvals, c(0.25, 0.75))
    qe <- quantile(lepvals, c(0.25, 0.75))
    slope <- diff(qp)/diff(qe)
    int <- qp[1] - slope * qe[1]
    if(plotline==1) {
      abline(0, 1, lwd = 1, col = colline)  
    } else if (plotline==2) {
      abline(int, slope, lwd = 1, col = colline)
    } else {
      abline(int, slope, lwd = 1, col = colline) # default
    } 
  }
  out <- list(pvals = pvals, epvals = epvals)
}

