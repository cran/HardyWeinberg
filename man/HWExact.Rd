\name{HWExact}
\alias{HWExact}
\title{ Exact test for Hardy-Weinberg equilibrium }
\description{
  \code{HWExact} performs an exact test for Hardy-Weinberg equilibrium
}
\usage{
HWExact(X, singleterms = FALSE, alternative = "two.sided", sorting =
"probability", verbose = FALSE)
}
\arguments{
  \item{X}{ vector of genotype counts AA,AB,BB }
  \item{singleterms}{ print probabilities and statistics for all
    possible samples (by default \code{FALSE})}
  \item{alternative}{ "two.sided" is the only option allowed so far}
  \item{sorting}{ sort samples by their probability (\code{sorting =
      "probability"}) or by the number of heterozygotes}
  \item{verbose}{print test results or not.}
}
\value{
  \item{D}{ deviation from independence }
  \item{pval}{ p-value of the exact test }
  \item{pofthesample}{ probability of the observed sample }
}
\references{Weir, B.S. (1996) Genetic data analysis II. Sinauer Associates, Massachusetts. See Chapter 3.}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{ \code{\link{HWLratio}} }
\examples{
x <- c(298,489,213)
names(x) <- c("MM","MN","NN")
HW.test <- HWExact(x)
}
\keyword{htest}

