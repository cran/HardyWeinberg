\name{HWChisqMat}
\alias{HWChisqMat}
\title{Matrix version of HWChisq}
\description{
   Function \code{HWChisqMat} executes the Chisquare test for HWE for each row in a matrix. This function is
   deprecated and it is better to use the faster \code{HWChisqStats} instead.
}
\usage{
HWChisqMat(X, ...)
}
\arguments{
  \item{X}{A n times 3 matrix of genotypic counts (AA,AB,BB)}
  \item{\dots}{extra arguments that are passed on to \code{HWChisqStats}}
}
\value{a vector with chi-square statistics or p-values}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWChisq}}, \code{\link{HWChisqStats}}  }
\examples{
X <- HWData(100,100)
colnames(X) <- c("MM","MN","NN")
Results <- HWChisqMat(X)
}
\keyword{htest}
