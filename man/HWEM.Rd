\name{HWEM}
\alias{HWEM}
\title{
Estimation of contributions of two populations to a sample of genotype frequencies with the EM algorithm.
}
\description{
Function \code{HWEM} estimates the relative contributions of two populations with different allele frequencies that contributed
to the genotype frequencies of sample which is a mixture of these two populations. 
}
\usage{
HWEM(x, p = NULL, G = NULL, delta.init = c(0.5, 0.5), itmax = 50,
     eps = 1e-06, verbose = FALSE)
}
\arguments{
  \item{x}{a 3x1 vector of genotype frequencies of the mixed sample, either counts or proportions.}
  \item{p}{a vector with two allele frequencies.}
  \item{G}{a 3x2 matrix with genotype frequencies as columns.}
  \item{delta.init}{initial values for the contributions of each populations}
  \item{itmax}{maximum number of iterations}
  \item{eps}{tolerance criterion for convergence}
  \item{verbose}{\code{verbose = TRUE} prints results, \code{verbose = FALSE} is silent.}
}
\details{
  \code{HWEM} employs the EM algorithm to estimate the contribution of two populations to the genotype frequencies of a bi-allelic
  marker, assumed to be typed for a mixed sample with individuals coming from two different populations with different allele
  frequencies. To estimate the contributions, both the mixed sample genotype frequencies (x) and either the genotype frequencies (G)
  or the allele frequencies (p) of the original populations must be specified. In case only allele frequencies are specified,
  Hardy-Weinberg proportions are assumed for the genotype frequencies. 
}
\value{a vector with two proportions, ordered according to the specified allele or genotype frequencies.}
\references{
Dempster, A.P., Laird, N.M. and Rubin, D.B. (1977) Maximum Likelihood from Incomplete Data via the EM Algorithm. Journal of the Royal
Statistical Society. Series B (Methodological) 39(1) pp. 1--38.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\examples{
# genotype frequencies population 1
g1 <- c(0.034, 0.330, 0.636)
# genotype frequencies population 2
g2 <- c(0.349, 0.493, 0.158)
# sample from the mixed population
x  <- c(0.270, 0.453,0.277)
#
# estimation based on genotype frequencies
#
G <- cbind(g1,g2)
contributions <- HWEM(x,G=G)
#
# estimation based on allele frequencies
#
p <- c(af(g1),af(g2))
contributions <- HWEM(x,p=p)
}
\keyword{misc}
