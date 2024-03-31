\name{CEUchr22}
\alias{CEUchr22}
\docType{data}
\title{
   Biallelic polymorphisms sampled from chromosome 22 of the CEU population of the 1000 Genomes project.
}
\description{
   Matrix \code{CEUchr22} contains 10000 single nucleotide polymorphisms sampled from chromosome 22 of sample of
   99 Utah residents with Northern and Western European ancestry.
}
\usage{data("CEUchr22")}
\format{Matrix}
\details{
   The polymorphisms are coded in (0,1,2) format representing the count of the reference allele.
}
\source{
   https://www.internationalgenome.org/
}
\references{
  The 1000 Genomes Project Consortium (2015) A global reference for human genetic variation. Nature 526(7571), pp. 68--74.
}
\examples{
data(CEUchr22)
str(CEUchr22)
}
\keyword{datasets}
