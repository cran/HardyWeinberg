\name{HWStrata}
\alias{HWStrata}
\title{
  Asymptotic test for HWE across strata for a single biallelic marker 
}
\description{
  Function \code{HWStrata} implements Olson's asymptotic test for HWE for a stratified sample of
  single biallelic polymorphism.
}
\usage{
HWStrata(X, verbose = TRUE)
}
\arguments{
  \item{X}{A three-column matrix of genotype counts (e.g. with columns AA, AB, BB)}
  \item{verbose}{print output if \code{verbose = TRUE}}
}
\details{
  See the references for the related homogeneity assumption.
}
\value{
  \item{T2}{The test statistic}
  \item{pval}{The p-value}
}

\references{
Olson J.M. (1993) Testing the Hardy-Weinberg law across strata. Annals of Human Genetics 57(4):291-295. 

Olson, J.M. and Foley, M. (1996) Testing for homogeneity of Hardy-Weinberg disequilibrium using data sampled from several
populations. Biometrics 52(3) pp. 971-979.
}

\author{ Jan Graffelman \email{jan.graffelman@upc.edu}} 

\seealso{
\code{\link{HWExactStats}}, \code{\link{HWChisqStats}} 
}

\examples{
#
# Test across strata
#
data("Glyoxalase")
Glyoxalase <- as.matrix(Glyoxalase)
HWStrata(Glyoxalase)
#
# Stratified exact testing, testing each sample
#
HWExactStats(Glyoxalase)
}
\keyword{htest}
