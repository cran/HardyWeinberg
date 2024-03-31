\name{af}
\alias{af}
\title{Function to compute allele frequencies}
\description{
  Function \code{af} computes the allele frequency for a vector with autosomal or X-chromosomal genotype counts or compositions.
}
\usage{
af(x)
}
\arguments{
  \item{x}{a vector or matrix with counts or compositions. \code{x} must have either three elements (autosomal) or five elements (X-chromosomal)}
}
\details{
Function \code{af} calculates the A allele frequency for the given order of autosomal genotype counts (AA,AB,BB) or the given order of
X chromosomal genotype counts (A,B,AA,AB,BB). The genotypes must be supplied this order. Use \code{maf} to calculate the minor allele
frequency and more flexible coding.
}
\value{
  a vector with allele frequencies
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{maf}}, \code{\link{afx}}}
\examples{
#
# A simulated random marker
#
X <- as.vector(rmultinom(1,100,c(0.5,0.4,0.1)))
X <- X/sum(X)
print(X)
print(af(X))
#
# MN blood group counts
#
x <- c(MM=298,MN=489,NN=213)
#
# Calculate frequency M allele
#
af(x)
#
# Calculate frequency N allele
#
af(rev(x))
1 - af(x)
#
# Calculated allele frequencies (p) for a matrix of genotype counts.
#
X <- HWData(nm=10,n=100)
X
p <- af(X)
p
}
\keyword{misc}

