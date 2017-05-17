\name{HWPosterior}
\alias{HWPosterior}
\title{
Calculation of posterior probabilities and Bayes factors for
Hardy-Weinberg tests at X-chromosomal variants.
}
\description{
Function \code{HWPosterior} calculates posterior probabilities and Bayes
factors for tests for Hardy-Weinberg equilibrium of X-chromosomal
variants. Four possible models are considered, and the posterior
probabilities and Bayes factors for each model are calculated. 
}
\usage{
HWPosterior(X, verbose = TRUE, prior.af = c(0.5,0.5), prior.gf = c(0.333,0.333,0.333))
}
\arguments{
  \item{X}{A vector of genotype counts. The order c(A,B,AA,AB,BB) is
    assumed. Differently ordered vectors can be supplied but then elements must
    be labeled by their genotype}
  \item{verbose}{prints results if \code{verbose = TRUE}}
  \item{prior.af}{Beta prior parameters for male and female allele frequencies}
  \item{prior.gf}{Dirichlet prior parameters for female genotype frequencies}
}
\details{
Default priors are used for genotype and allele frequencies.
}
\value{
  A matrix with posterior probabilities and Bayes factors.
}
\references{
  Puig, X., Ginebra, J. and Graffelman, J. (2017) A Bayesian test for
  Hardy-Weinberg equilibrium of bi-allelic X-chromosomal markers. To
  appear in Heredity.
}
\author{
  Xavi Puig \email{xavier.puig@upc.edu} and
  Jan Graffelman \email{jan.graffelman@upc.edu}
}
\seealso{ \code{\link{HWChisq}}, \code{\link{HWExact}}, \code{\link{HWExactStats}} }
\examples{
x <- c(A=43,B=13,AA=26,AB=19,BB=3)
out <- HWPosterior(x,verbose=TRUE)
}
\keyword{htest}

