\name{HWChisq}
\Rdversion{1.1}
\alias{HWChisq}
\title{Chi square tests for Hardy Weinberg equilibrium}
\description{
  \code{HWChisq} performs the chi-square test for Hardy Weinberg
  equilibrium with or without continuity correction.
}
\usage{
HWChisq(X, cc = 0.5, alpha = 0.05, verbose = FALSE)
}
\arguments{
  \item{X}{\code{X} a vector containg the genotypic counts (AA,AB,BB).}
  \item{cc}{\code{cc} the continuity correction parameter (default
    \code{cc = 0.5}).}
  \item{alpha}{ significance level (0.05 by default). }  
  \item{verbose}{ \code{verbose} = 1 prints results, \code{verbose} = 0 is silent.}
}
\details{
   \code{HWChisq} does a chi-square test for Hardy-Weinberg equilibrium,
   and by default applies a continuity correction. For extreme allele
   frequencies, the continuity correction can lead to excessive type 1
   error rates, and is better turned off in that case. The continuity
   correction can be turned off by specifying \code{cc=0}.
}
\value{
\code{HWChisq} returns a list with the components:
  \item{chisq }{value of the chi-square statistic. NA is returned if the marker is monomorphic.}
  \item{pval }{p-value of the chi-square test for Hardy-Weinberg equilibrium.}
  \item{D }{Half the deviation from Hardy-Weinberg equilibrium for the AB genotype.}
  \item{p }{the allele frequency of A.}
  \item{f }{the inbreeding coefficient.}
}
\references{Weir, B.S. (1996) Genetic data analysis II. Sinauer
  Associates, Massachusetts. See Chapter3.}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWLratio}} }
\examples{
x <- c(MM=298,MN=489,NN=213)
HW.test <- HWChisq(x,verbose=TRUE)
}
\keyword{htest}
