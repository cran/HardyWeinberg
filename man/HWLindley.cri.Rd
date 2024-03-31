\name{HWLindley.cri}
\alias{HWLindley.cri}
\title{
  Calculate a credible interval for Lindley's alpha for HWE,
}
\description{
  Function \code{HWLindley.cri} calculates a Bayesian credible interval using Lindley's
  posterior density for equilibrium paramater alpha.
}
\usage{
HWLindley.cri(x, verbose = TRUE, limits = c(0.025, 0.975))
}
\arguments{
  \item{x}{a vector of three genotype counts in order (AA,AB,BB).}
  \item{verbose}{print output (\code{verbose=TRUE}) or be silent}
  \item{limits}{upper and lower probability limits of the interval}
}
\details{
  The limits are found by numerical integration over Lindley's density.
}
\value{a vector with the lower and upper limit of the credible interval}
\references{
  Lindley, D.V. (1988) Statistical Inference Concerning Hardy-Weinberg Equilibrium. In: 
  Bernardo, J.M., DeGroot, M.H., Lindley, D.V. and Smith, A.F.M. Bayesian Statistics, 3,
  pp. 307-326. Oxford University Press.

  Graffelman, J. (2020) Statistical tests for the Hardy-Weinberg equilibrium.
  Wiley StatsRef: Statistics Reference Online \doi{10.1002/9781118445112.stat08274}.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{
  \code{\link{HWLindley}}, \code{\link{HWPosterior}},
}
\examples{
#
# MN blood group data
#
x <- c(MM=298,MN=489,NN=213)
#
# credible interval of 95%
#
HWLindley.cri(x)
#
# credible interval of 90%
#
HWLindley.cri(x,limits=c(0.05,0.95))
}
\keyword{htest}
