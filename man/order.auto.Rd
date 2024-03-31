\name{order.auto}
\alias{order.auto}
\title{Reordering of autosomal genotype counts}
\description{
  Function \code{order.auto} reorders a named vector of three genotype
  counts, such that the sequence (minor homozygote, heterozygote,
  major heterzygote) is establised.
}
\usage{order.auto(X)}
\arguments{
  \item{X}{a named vector of genotype counts (e.g. c(AA=25,AB=50,BB=25))}
}
\value{a vector}
\references{
Graffelman, J. (2015) Exploring Diallelic Genetic Markers: The HardyWeinberg Package.
\emph{Journal of Statistical Software} 64(3): 1-23. \doi{10.18637/jss.v064.i03}.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{
  \code{\link{order.x}}
}
\examples{
x <- c(MN=489,MM=298,NN=213)
print(x)
y <- order.auto(x)
print(y)
}
\keyword{manip}

