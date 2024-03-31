\name{HWExactMat}
\alias{HWExactMat}
\title{Matrix version of HWExact}
\description{
   Function \code{HWExactMat} is deprecated; use \code{HWExactStats} instead. 
}
\usage{
HWExactMat(X, ...)
}
\arguments{
  \item{X}{A n times 3 matrix of genotypic counts (AA,AB,BB)}
  \item{\dots}{extra arguments that are passed on to \code{HWExact}}
}
\value{
  \item{pvalvec}{Vector with the p-values of each test}
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWExact}} \code{\link{HWExactStats}} }
\examples{
X <- HWData(100,100)
colnames(X) <- c("MM","MN","NN")
Results <- HWExactMat(X)
Output <- cbind(X,Results)
print(Output)
}
\keyword{htest}
