\name{HWExact}
\alias{HWExact}
\title{ Exact test for Hardy-Weinberg equilibrium }
\description{
  \code{HWExact} performs an exact test for Hardy-Weinberg equilibrium
}
\usage{
HWExact(X, alternative = "two.sided", pvaluetype = "selome", verbose = FALSE)
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
equally likely as the current sample. When set to \code{midp}, the
p-value is computed as half the probability of the current sample +
the probabilities of all samples that are more extreme.}
   \item{verbose}{ print results or not.}
}
\details{
   \code{HWExact} uses the recursion equations described by Wigginton
   et. al.

   For large samples, \code{HWExact} may give the error message:
   "evaluation nested too deeply: infinite recursion". This can usually
   be resolved by increasing R's limit on nested expressions with
   \code{options(expressions=10000)} or a higher limit. With higher
   limits, the error message "protect(): protection stack overflow" can
   occur. This error can usually be resolved by increasing R's
   protection stack with the command line option
   \code{--max-ppsize 100000} or higer values. However, with such large
   samples the exact test will give virtually the same result as a
   chi-square test, and it may be easier to use \code{HWChisq} in these
   circumstances.
}
\value{
  \item{pval}{ p-value of the exact test }
  \item{prob}{ probabilities of all possible samples with the same
   sample size and minor allele count}
  \item{pofthesample}{ probability of the observed sample }
}
\references{Weir, B.S. (1996) Genetic data analysis II. Sinauer
  Associates, Massachusetts. See Chapter3.
  
  Wigginton, J.E., Cutler, D.J. and Abecasis, G.R. (2005) A note on exact
  tests of Hardy-Weinberg equilibrium, American Journal of Human
  Genetics (76) pp. 887-893.}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{ \code{\link{HWLratio}}, \code{\link{HWChisq}} }
\examples{
x <- c(298,489,213)
names(x) <- c("MM","MN","NN")
HW.test <- HWExact(x)
print(HW.test)
}
\keyword{htest}
