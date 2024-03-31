\name{order.x}
\alias{order.x}
\title{Reordering of X-chromosomal genotype counts}

\description{
  Function \code{order.x} reorders a named vector of three genotype
  counts, such that the sequence (minor hemizygote, major hemizygote,
  minor homozygote, heterozygote, major heterzygote) is establised.
}
\usage{order.x(X)}
\arguments{
  \item{X}{a named vector of genotype counts (e.g. c(A=10,B=10,AA=25,AB=50,BB=25))}
}
\value{a vector}
\references{
Graffelman, J. (2015) Exploring Diallelic Genetic Markers: The HardyWeinberg Package.
\emph{Journal of Statistical Software} 64(3): 1-23. \doi{10.18637/jss.v064.i03}.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{
  \code{\link{order.auto}}
}
\examples{
x <-  c(A=392, B=212, AA=275, AB=296, BB=80)
print(x)
y <- order.x(x)
print(y)
}
\keyword{manip}

