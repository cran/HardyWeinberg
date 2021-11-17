\name{AllelesToTriangular}
\alias{AllelesToTriangular}
\title{
  Calculate triangular genotype matrix for vector(s) of alleles.
}
\description{
  \code{AllelesToTriangular} constructs a lower triangular matrix of genotype counts from one or two vectors of alleles. It is particularly useful to create genotype counts for microsatellite data (STRs). 
}
\usage{
AllelesToTriangular(A1, A2 = NULL, given=NULL) 
}
\arguments{
  \item{A1}{The first allele of each individual, or a vector with all alleles, two consecutive ones for each individual.}
  \item{A2}{The second allele of each individual (optional).}
  \item{given}{A vector of known alleles (optional). This argument can be used to specify alleles that may not exist in the data.}
}
\details{
  If the data is a single column vector with two succesive alleles for each individual, then specify \code{A1} only. If data consists of two columns, each holding one allele of each individual, then specify \code{A1} and \code{A2}. Typical STR data that comes in the format of two repeat lengths for a set of individuals can be transformed into a lower triangular matrix with genotype counts. See the examples below.
}
\value{
  A lower triangular matrix with genotype counts.
}
\references{
Graffelman, J. (2015) Exploring Diallelic Genetic Markers: The HardyWeinberg Package.
\emph{Journal of Statistical Software} 64(3): 1-23. \doi{10.18637/jss.v064.i03}.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{ \code{\link{toTriangular}} }
\examples{
data(NistSTRs)
A1 <- NistSTRs[,1]
A2 <- NistSTRs[,2]
GM <- AllelesToTriangular(A1,A2)
print(GM)
}
\keyword{misc}

