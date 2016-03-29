\name{HWExactStats}
\alias{HWExactStats}
\title{
Computation of Exact p-values for Hardy-Weinberg equilibrium for sets of SNPs
}
\description{
  \code{HWExactStats} is a function for the computation of Exact p-values
  for a large set of bi-allelic markers (typically SNPs). 
}
\usage{
HWExactStats(X, x.linked = FALSE, ...)
}
\arguments{
  \item{X}{A matrix with genotype counts, one row per marker. \code{X}
    should have 5 columns for an X-chromosomal data set and 3 columns
    for an autosomal data set.}
  \item{x.linked}{Logical indicating whether the markers are autosomal (\code{x.linked=FALSE})
    or X-chromosomal (\code{x.linked=TRUE}).}
  \item{\dots}{extra arguments that are passed on to \code{HWExact}}
}
\details{
Matrix \code{X} should strictly comply with the following format. For
an autosomal dataset is should contain the 3 genotype counts in order
(AA,AB,BB). For an X-chromosomal dataset it should contain the 5
genotype counts in order (A,B,AA,AB,BB) where A and B are the male
counts and AA, AB and BB the female counts.
}
\value{A vector of p-values}
\references{
  Graffelman, J. and Weir, B.S. (2016) Testing for Hardy-Weinberg
  equilibrium at 2 bi-allelic genetic markers on the X chromosome.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{HWExact}}  }
\examples{
#
# Autosomal example
#
set.seed(123)
X <- HWData(1000,100)
monom <- (X[,2]==0 & X[,1]==0) | (X[,2]==0 & X[,3]==0)
X <- X[!monom,] # exclude monomorphics
Exact.pvalues <- HWExactStats(X,x.linked=FALSE)

#
# X-chromosomal example
#
X <- HWData(1000,100,n.males=50,nA=75,x.linked=TRUE)
Exact.pvalues <- HWExactStats(X,x.linked=TRUE)
}
\keyword{htest}
