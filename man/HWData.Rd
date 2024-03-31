\name{HWData}
\alias{HWData}
\title{Generate genetic marker data in or out of Hardy-Weinberg Equilibrium}
\description{
   HWData generates samples of genotypic counts under various schemes. It mainly
   uses sampling from the multinomial distribution for given or random allele frequencies,
   either assuming Hardy-Weinberg proportions or a specified degree of inbreeding. Sampling
   can also be performed conditional on the allele frequency. The same procedures are
   also available for X linked markers.
}
\usage{
HWData(nm = 100, n = rep(100, nm), f = rep(0, nm), p = NULL, conditional = FALSE,
exactequilibrium = FALSE, x.linked = FALSE, nA = NULL, n.males = round(0.5 * n),
shape1 = 1, shape2 = 1, counts = TRUE)
}
\arguments{
  \item{nm}{The number of bi-allelic markers.}
  \item{n}{The sample sizes.}    
  \item{f}{The inbreeding coefficients (only for autosomal markers)}
  \item{p}{a vector of allele frequencies}
  \item{conditional}{if \code{TRUE} The Levene-Haldane distribution will be used for sampling autosomal markers,
   the Graffelman-Weir distribution for sampling X chromosomal markers. If \code{FALSE} a multinomial distribution is used}
  \item{exactequilibrium}{generates data in exact HWE if set to
    \code{TRUE}}
  \item{x.linked}{Simulated autosomal markers (\code{x.linked=FALSE},
    the default) or X-chromosomal markers (\code{x.linked=TRUE})}
  \item{nA}{A vector of minor allele counts, one for each marker. If not
    specified, it will be calculated from \code{p}}
  \item{n.males}{The number of males (only relevant if \code{x.linked = TRUE})}
  \item{shape1}{First shape parameter of the beta distribution used to generate allele frequencies}
  \item{shape2}{Second shape parameter of the beta distribution used to generate allele frequencies}
  \item{counts}{If \code{counts=TRUE} the output are genotype counts, if \code{counts=TRUE} relative genotype frequencies}  
}
\details{
  Option \code{pfixed} is deprecated and replaced by \code{conditional}
  Option \code{pdist} is deprecated and replaced by parameters \code{shape1} and \code{shape2}

  \code{HWData} returns a matrix of genotype counts, \code{nm} by 3 for
  autsomal markers or \code{nm} by 5 for X-chromosomal markers.

  If the inbreeding coefficient is specified (\code{f}) it will only
  take effect for autosomal markers (\code{x.linked=FALSE}) and
  multinomial sampling (\code{conditional=FALSE}).
}
\value{
  \item{X}{A matrix containing the genotype counts.}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{HWTernaryPlot}}}
\examples{
#
# Generate 100 SNPs with uniform allele frequency under the equilibrium assumption.
#
out <- HWData(nm=100,n=100)
#
# Generate genotype frequencies of 100 SNPs with uniform allele frequency assuming exact equilbrium.
#
X <- HWData(nm=100, exactequilibrium = TRUE, counts = FALSE)
#
# Generate 100 SNPs (as counts), all having an expected A allele frequency of 0.50
#
X <- HWData(nm=100,p=0.5)
#
#  Generate 100 SNPs, 50 with A allele frequency 0.25 and 50 with A allele frequency 0.75,
#  assuming fixed allele frequencies. 
#
X <- HWData(nm=100,p=rep(c(0.25,0.75),50), conditional = TRUE)
#
#  Generate 100 SNPs with a skewed (beta) distribution of the allele frequency,
#  rich in variants with a low minor allele frequency.
#
X <- HWData(nm=100,shape1=1,shape2=10)
#
# Generate 100 X chromosomal SNPs with uniformly distributed allele frequency.
#
X <- HWData(nm=100, x.linked = TRUE)
}
\keyword{datagen}
