\name{is.mono}
\alias{is.mono}
\title{Detects monomorphic variants}
\description{
Function \code{is.mono} tests if a bi-allelic variant is monomorphic or not}
\usage{
is.mono(x)
}
\arguments{
  \item{x}{a vector of three genotype counts (AA,AB,BB), or a three-column matrix with
    genotype counts (variants in rows, genotype counts (AA,AB,BB) in three columns)}
}
\value{A logical or vector of logicals}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
#
# a polymorphic marker
#
x <- c(AA=10,AB=20,BB=10)
print(is.mono(x))
#
# a monomorphic marker
#
x <- c(AT=0,AA=100,TT=0)
print(is.mono(x))
#
# a marker with only heterozygotes
#
x <- c(AT=100,AA=0,TT=0)
print(is.mono(x))
#
# a matrix with low maf markers
#
set.seed(123)
X <- HWData(100,50,shape1=1,shape2=20)
number.monomorphics <- sum(is.mono(X))
print(number.monomorphics)
}
\keyword{misc}
