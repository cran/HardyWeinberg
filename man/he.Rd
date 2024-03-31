\name{he}
\alias{he}
\title{
  Calculate expected heterozygosity (He)
}
\description{
  Function \code{he} calculates the expected heterozygosity for a biallelic genetic variant.
}
\usage{
he(x, bias.correct = TRUE)
}
\arguments{
  \item{x}{a vector or matrix with genotype counts.}
  \item{bias.correct}{if \code{bias.correct = TRUE} a correction for bias will be applied.}
}
\details{
  \code{x} can be a vector of genotype counts for a single marker (AA,AB,BB) or a three-column matrix of
  genotype counts for multiple markers.
}
\value{
 a vector of expected heterozygosities.
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
  \code{\link{af}}, \code{\link{maf}}.
}
\examples{
#
# He for a single marker
#
x <- c(MM=298,MN=489,NN=213)
he(x)
#
# He for a matrix of rmarkers
#
set.seed(123)
X <- HWData(10,100)
he(X)
he(X,bias.correct = FALSE)
}
\keyword{univar}
