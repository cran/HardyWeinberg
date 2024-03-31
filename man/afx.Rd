\name{afx}
\alias{afx}
\title{Function to compute X-chromosomal allele frequencies}
\description{
  Function \code{afx} computes the allele frequency for a vector or matrix of X-chromosomal genotype counts or frequencies.
}
\usage{
afx(x)
}
\arguments{
  \item{x}{a vector or matrix with X-chromosomal counts or compositions. \code{x} must have five elements or columns.}
}
\details{
Function \code{afx} calculates the A allele frequency for the given order of X-chromosomal genotype counts (A,B,AA,AB,BB). The genotypes must be supplied this order. Use \code{maf} to calculate the minor allele frequency and more flexible coding.
}
\value{
  a vector with allele frequencies
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{maf}}, \code{\link{af}}}
\examples{
#
# A simulated random X-chromosomal marker
#
X <- as.vector(rmultinom(1,100,c(0.3,0.2,0.1,0.15,0.25)))
X <- X/sum(X)
print(X)
print(afx(X))
#
# 
#
rs5968922 <-  c(A=392, B=212, AA=275, AB=296, BB=80)
afx(rs5968922)
#
#
#
y <- c(C=337, G=129, CC=277, CG=209, GG=48) 
afx(y)
#
# Calculated allele frequencies (p) for a matrix of X-chromosomal genotype counts.
#
X <- HWData(nm=10,n=1000,x.linked=TRUE)
X
p <- afx(X)
p
}
\keyword{misc}

