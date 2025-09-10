\name{HWBalance}
\alias{HWBalance}
\title{
Compositional Statistic and Test for Hardy-Weinberg equilibrium
}
\description{
Function \code{HWBalance} can test bi and multiallelic genetic polymorphisms for Hardy-Weinberg equilibrium using the ratio of the
geometric means of heterozygote and homozygote genotype frequencies.
}
\usage{
HWBalance(x, verbose = TRUE, alpha = 0.05, alternative = "two.sided")
}
\arguments{
  \item{x}{A 3x1 vector of genotype counts or a lower triangular matrix of genotype counts.}
  \item{verbose}{\code{verbose = TRUE} prints results, \code{verbose = FALSE} is silent.}
  \item{alpha}{The significance level.}
  \item{alternative}{\code{alternative="two.sided"} tests for equilibrium; \code{alternative="greater"} tests against the alternative of heterozygote excess; \code{alternative="smaller"} the alternative of a lack of heterozygotes.}
}
\details{
  The test is based on an isometric logratio, the ratio of the geometric means of heterozygotes and homozygotes. This ratio has asymptotically a normal distribution. Note that this test does not allow for genotype counts that are zero.    
}
\value{
   \item{Z1}{The value of the heterozygote-homozygote balance.}
   \item{VZ1}{The sample variance of the heterozygote-homozygote balance.}
   \item{Z}{The Z statistic (to be referred to a standard normal distribution).}
   \item{chisq}{The Chi-square statistic.}
   \item{pval}{The p-value.}
}
\references{
Graffelman, J. (2025) A logratio approach to the analysis of autosomal genotype frequencies across multiple samples. Preprint at www.bioRxiv.org. \doi{10.1101/2025.01.17.633675}
}

\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }

\seealso{ \code{\link{HWLratio}}, \code{\link{HWChisqStats}}, \code{\link{HWExact}} }

\examples{
#
# Test for an autosomal biallelic blood group marker
#
x <- c(MM=298,MN=489,NN=213)
HW.test <- HWBalance(x,verbose=TRUE)
#
#
#
}

\keyword{htest}

