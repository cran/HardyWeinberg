\name{HWExact}
\alias{HWExact}
\title{ Exact test for Hardy-Weinberg equilibrium }
\description{
  \code{HWExact} performs an exact test for Hardy-Weinberg equilibrium
}
\usage{
HWExact(X, alternative = "two.sided", pvaluetype = "dost", verbose = FALSE)
}
\arguments{
  \item{X}{ vector with the genotype counts AA, AB, BB }
  \item{alternative}{\code{two.sided} (default) will perform a two-sided
    test where both an excess and a dearth of heterozygotes count as
    evidence against HWE. \code{less} is a one-sided test where only
    dearth of heterozygotes counts a evidence against HWE,
    \code{greater} is a one-sided test where only excess of
    heterozygotes counts as evidence against HWE.}
  \item{pvaluetype}{if \code{pvaluetype} is set to \code{dost}
then the p-value of a two-sided test is computed as twice the tail area
of a one-sided test. When set to \code{selome}, the p-value
is computed as the sum of the probabilities of all samples less or
equally likely as the current sample.}
   \item{verbose}{ print results or not.}
}
\details{
   \code{HWExact} use the recursion equations described by Wigginton
   et. al.
}
\value{
  \item{pval}{ p-value of the exact test }
  \item{pofthesample}{ probability of the observed sample }
}
\references{Weir, B.S. (1996) Genetic data analysis II. Sinauer
  Associates, Massachusetts. See Chapter3.
  
  Wigginton, J.E., Cutler, D.J. and Abecasis, G.R. (2005) A note on exact
  tests of Hardy-Weinberg equilibrium, American Journal of Human
  Genetics (76) pp. 887-893.}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\note{
   \code{HWExact} is designed for a fast analysis of a large set of
   markers. \code{HWExact} is slower, but provides more detailed
   statistics for the sample under study.
}
\seealso{ \code{\link{HWLratio}}, \code{\link{HWChisq}} }
\examples{
x <- c(298,489,213)
names(x) <- c("MM","MN","NN")
HW.test <- HWExact(x)
print(HW.test)
}
\keyword{htest}
