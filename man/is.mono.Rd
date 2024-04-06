\name{is.mono}
\alias{is.mono}
\title{Detects autosomal and X-chromosomal monomorphic variants}
\description{
Function \code{is.mono} tests if a bi-allelic variant is monomorphic or not}
\usage{
is.mono(x)
}
\arguments{
  \item{x}{a vector of three or five genotype counts ((AA,AB,BB) or (A,B,AA,AB,BB)), or a
           three-column or five-column matrix with genotype counts (variants in rows, columns)}
}
\details{
   \code{is.mono} assumes autosomal variants are coded in three-element vectors or
   three-column matrices, whereas X-chromosomal variants are coded in five-element
   vectors or five-column matrices.
}
\value{A logical or vector of logicals}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
#
# a polymorphic autosomal marker
#
x <- c(AA=10,AB=20,BB=10)
print(is.mono(x))
#
# a monomorphic autosomal marker
#
x <- c(AT=0,AA=100,TT=0)
print(is.mono(x))
#
# an autosomal marker with only heterozygotes
#
x <- c(AT=100,AA=0,TT=0)
print(is.mono(x))
#
# a matrix with low maf autosomal markers
#
set.seed(123)
X <- HWData(100,50,shape1=1,shape2=20)
number.monomorphics <- sum(is.mono(X))
print(number.monomorphics)
#
# a polymorphic X chromosomal marker
#
x <- c(G=24,C=26,GG=12,CC=13,GC=25)
print(is.mono(x))
#
# another polymorphic X chromosomal marker
#
x <- c(G=24,C=1,GG=25,CC=0,GC=0)
is.mono(x)
#
# a monomorphic X chromosomal marker
#
x <- c(G=24,C=0,GG=12,CC=0,GC=0)
is.mono(x)
#
# a matrix with low maf X-chromosomal markers
#
set.seed(123)
Y <- HWData(100,50,shape=1,shape2=20,x.linked = TRUE)
number.monomorphics <- sum(is.mono(Y))
print(number.monomorphics)
}
\keyword{misc}
