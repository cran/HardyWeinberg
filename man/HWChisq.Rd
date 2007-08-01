\name{HWChisq}
\alias{HWChisq}
\title{Chi square tests for Hardy Weinberg equilibrium}
\description{
  \code{HWChisq} performs the chi-square test for Hardy Weinberg
  equilibrium with or without continuity correction.
}
\usage{
HWChisq(X, cc = 0, verbose = FALSE)
}
\arguments{
  \item{X}{\code{X} a vector containg the genotypic counts (AA,AB,BB).}
  \item{cc}{\code{cc} the continuity correction parameter (default \code{cc = 0}).}
  \item{verbose}{ \code{verbose} = 1 prints results, \code{verbose} = 0 is silent.}
}
\value{
\code{HWChisq} returns a list with the components:
  \item{chisq }{value of the chi-square statistic. NA is returned if the marker is monomorphic.}
  \item{pval }{p-value of the chi-square test for Hardy-Weinberg equilibrium.}
  \item{D }{Half the deviation from Hardy-Weinberg equilibrium for the AB genotype.}
  \item{p }{allele frequency of A.}
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWLratio}} }
\examples{
x <- c(298,489,213)
names(x) <- c("MM","MN","NN")
HW.test <- HWChisq(x,verbose=TRUE)
}
\keyword{htest}
