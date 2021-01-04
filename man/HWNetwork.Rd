\name{HWNetwork}
\alias{HWNetwork}
\title{
  Autosomal and X-chromosomal exact tests for HWE via a Network algorithm 
}
\description{
  Program \code{HWNetwork} implements a network algorithm for efficient calculation of exact test p-values in HWE tests with multiple alleles.
}
\usage{
HWNetwork(a1, a2, ma = NULL, fe = NULL, gender = NULL, verbose = TRUE)
}
\arguments{
  \item{a1}{the first allele (expressed as a number)}
  \item{a2}{the second allele (expressed as a number; NA if the variant is X chromosomal)}
  \item{ma}{alternative format: vector of male X chromosomal allele counts.}
  \item{fe}{triangular matrix of female genotype counts}
  \item{gender}{gender of the individual (1=male; 2=female)}
  \item{verbose}{be silent (\code{verbose=FALSE}) or informative (\code{verbose=TRUE})}
}
\details{
  Function \code{HWNetwork} accepts data in two formats. Original genotype data (e.g. repeat numbers of microsatellites) can be supplied, or  the data can be supplied in summarized form as a male genotype count vector and a female genotype count matrix. If one of the two male alleles is missing (NA) the variant will taken to be X-chromosomal. If all males have two alleles, the variant will taken to be autosomal.
}
\value{
  the exact p-value of the test.
}
\references{
Aoki, S. (2003) Network algorithm for the Exact Test of Hardy-{W}einberg Proportion for Multiple Alleles. \emph{Biometrical Journal} 45(4),
pp. 471-490.

Engels, W. R. (2009) Exact Tests for Hardy-Weinberg Proportions. \emph{Genetics} 183, pp. 1431-1441.

Graffelman, J. (2015) Exploring Diallelic Genetic Markers: The HardyWeinberg Package.
\emph{Journal of Statistical Software} 64(3): 1-23. \url{https://www.jstatsoft.org/v64/i03/}.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu}} 
\seealso{
 \code{\link{HWExact},\link{HWExactStats}}
}
\examples{
#
# From vectors with counts of genotypes
#
data(TSIXTriAllelics)
ma <- as.matrix(TSIXTriAllelics[1,2:4])
names(ma) <- c("A","B","C")

fe <- TSIXTriAllelics[1,5:10]
names(fe) <- c("AA","AB","AC","BB","BC","CC")

fe <- toTriangular(fe)

HWNetwork(ma=ma,fe=fe)
}
\keyword{htest}
