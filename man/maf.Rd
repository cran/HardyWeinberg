\name{maf}
\alias{maf}
\title{Function to compute minor allele frequencies}
\description{
  Function \code{maf} computes the minor allele frequency for
  a matrix or vector of genotype counts.
}
\usage{
maf(x, option = 1, verbose = FALSE) 
}
\arguments{
  \item{x}{a vector or matrix of with genotype counts}
  \item{option}{determines output that is returned. \code{option=1} calculates the minor allele frequency; \code{option=2}
  returns all allele frequencies; \code{option=3} returns the allele counts.}
  \item{verbose}{be silent (\code{verbose=FALSE}) or not.}
}
\value{
  a vector of minor allele frequencies or minor allele counts.
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{
 \code{\link{mac}}, \code{\link{MakeCounts}}
}
\examples{
#
# MAF of a single random marker
#
set.seed(123)
X <- as.vector(rmultinom(1,100,c(0.5,0.4,0.1)))
names(X) <- c("AA","AB","BB")
print(X)
print(maf(X))
#
# MAF of MN bloodgroup counts
#
x <- c(MM=298,MN=489,NN=213)
maf(x)
#
# Both allele frqeuencies, ordered from minor to major
#
maf(x,2)
#
# allele counts of MN bloodgroup counts, order from minor to major
#
maf(x,3)
#
# MAF of single triallelic marker in triangular format
#
x <- c(AA=20,AB=52,AC=34,BB=17,BC=51,CC=26)
print(x)
GT <- toTriangular(x)
print(GT)
maf(GT)
#
# extract all allele frequencies
#
maf(GT,option=2)
#
# extract all allele counts
#
maf(GT,option=3)
#
# Calculate the MAF for 10 random SNPs under HWE
#
set.seed(123)
Z <- HWData(nm=10)
print(Z)
#
# vector with minor allele frequencies, one per marker
#
maf(Z)
#
# Matrix with minor and major allele frequencies, one row per marker 
#
maf(Z,2)
#
# Matrix with minor and major allele count, one row per marker 
#
maf(Z,3)
}
\keyword{misc}

