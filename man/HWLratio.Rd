\name{HWLratio}
\alias{HWLratio}
\title{Likelihood ratio test for Hardy Weinberg equilibrium}
\description{
  \code{HWLratio} performs the Likelihood ratio test for Hardy Weinberg equilibrium.
}
\usage{
HWLratio(X, verbose = FALSE)
}
\arguments{
  \item{X}{\code{X} a vector containing the genotypic counts (AA,AB,BB).}
  \item{verbose}{\code{verbose} = 1 prints results, \code{verbose} = 0 is silent.}
}
\value{
\code{HWLratio} returns a list with the components:
  \item{Lambda }{the likelihood ratio}
  \item{G2 }{-2*log(Lambda)}
  \item{pval}{the p-value}
}
\references{Weir, B.S. (1996) Genetic data analysis II. Sinauer Associates, Massachusetts. See Chapter 3.}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWChisq}} }
\examples{
x <- c(298,489,213)
names(x) <- c("MM","MN","NN")
HW.test <- HWLratio(x,verbose=TRUE)
}
\keyword{htest}
