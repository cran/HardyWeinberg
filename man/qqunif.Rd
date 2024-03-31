\name{qqunif}
\alias{qqunif}
\title{
Q-Q plot for a uniform distribution.
}
\description{
   \code{qqunif} makes a Q-Q plot against a uniform distribution for the supplied data vector.
}
\usage{
qqunif(x, logplot = FALSE, lbs = 1:length(x), texton = FALSE, 
          xylim = NULL, main = "Q-Q plot for a uniform distribution",
          plotline = 0, xlab = "Expected p-value", 
          ylab = "Observed p-value", colvec=rep("black",length(x)),
	  colline = "black", ...)
}
\arguments{
  \item{x}{The data vector}
  \item{logplot}{If \code{logplot} is set to true, then the log10 of the
  p-values will be used in the plot. If not, untransformed p-values will
  be used.}
  \item{lbs}{A vector of labels for the points in the Q-Q plot}
  \item{texton}{Logical indicating whether labels should be plotted or not}
  \item{xylim}{Shared upper limit x and y axis in the plot}
  \item{main}{Title for the plot}
  \item{plotline}{Setting \code{plotline=0} plots a (0,1) line; \code{plotline=1} plots a robust fit,}
  \item{xlab}{label for the x axis}
  \item{ylab}{label for the y axis}
  \item{colvec}{vector with colours for the points in the QQ plot}
  \item{colline}{colour for the line in the plot}
  \item{\dots}{Any additional arguments for the \code{plot} instruction}
}
\value{
   \item{pvals}{observed probabilities}
   \item{epvals}{expected probabilities}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{qqnorm}}}
\examples{
   x <- runif(1000)
   z <- qqunif(x)
}
\keyword{aplot}

