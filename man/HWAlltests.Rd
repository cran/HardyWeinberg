\name{HWAlltests}
\alias{HWAlltests}
\title{
Perform all tests for Hardy-Weinberg equilibrium
}
\description{
  \code{HWAlltests} performs all classical frequentists tests for
  Hardy-Weinberg equilibrium and lists their p-values.
}
\usage{
HWAlltests(x, verbose = FALSE, include.permutation.test=FALSE)
}
\arguments{
  \item{x}{a vector with a set of genotype counts (AA, AB, BB)}
  \item{verbose}{print output if set to TRUE}
  \item{include.permutation.test}{turns on the permutation test if set to TRUE}
}
\details{
By default the permutation test is not performed in order to reduce
computing time. 
   }
\value{
A dataframe with test statistics and p-values. 
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWLratio}}, \code{\link{HWChisq}}, \code{\link{HWExact}}  }
\examples{
x <- c(298,489,213)
names(x) <- c("MM","MN","NN")
HWAlltests(x,verbose=TRUE)
}
\keyword{htest}
