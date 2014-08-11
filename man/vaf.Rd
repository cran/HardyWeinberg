\name{vaf}
\alias{vaf}
\title{
Computes the sample variance of the minor allele frequencies
}
\description{
Function \code{vaf} computes the sample variance of the minor allele frequencies
of a single sample or a matrix of samples.
}
\usage{
vaf(X, hw = FALSE, minor = TRUE)
}
\arguments{
  \item{X}{vector or matrix with genotype counts (AA,AB,BB)}
  \item{hw}{assume Hardy-Weinberg proportions (\code{hw=TRUE}) or not (\code{hw=FALSE})}
  \item{minor}{compute the variance for the minor allele \code{minor=TRUE} or major allele \code{major=TRUE}}
}
\value{
  a numeric vector of variances.
}
\references{Weir, B.S. (1996) Genetic data analysis II. Sinauer Associates, Massachusetts. See Chapter 2.}
\author{Jan Graffelman \email{jan.graffelman@upc.edu}}
\seealso{\code{\link{af}},\code{\link{maf}}}
\examples{
x <- c(MM=298,MN=489,NN=213)
pA <- af(x)
vA <- vaf(x)
cat("allele frequency:",pA,"\n")
cat("sample variance allele frequency:",vA,"\n")
}
\keyword{misc}