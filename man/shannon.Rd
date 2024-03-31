\name{shannon}
\alias{shannon}
\title{
Shannon index
}
\description{
Function \code{shannon} calculates the Shannon index and its variance for a vector of counts.
}
\usage{
shannon(x)
}
\arguments{
  \item{x}{a vector of counts}
}
\value{
  \item{Hp}{the sample Shannon index}
  \item{VHp}{the sample variance of the Shannon index}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\examples{
#
# Shannon index for allele frequencies of a biallelic MN blood group polymorphism
#
x <- c(MM=298,MN=489,NN=213)
p <- af(x)
shannon(c(p,1-p))$Hp
#
# Shannon index for the allele frequencies of an STR
#
data("NistSTRs")
AlleleTable <- table(c(NistSTRs[,1],NistSTRs[,2]))
AlleleFreq  <- AlleleTable/sum(AlleleTable)
shannon(AlleleFreq)$Hp
}
\keyword{misc}% use one of  RShowDoc("KEYWORDS")

