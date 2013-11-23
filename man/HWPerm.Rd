\name{HWPerm}
\alias{HWPerm}
\title{
Permutation test for Hardy-Weinberg equilibrium
}
\description{
Function \code{HWPerm} does a permutation test for Hardy-Weinberg
equilibrium using a user-supplied test statistic.
}
\usage{
HWPerm(x, nperm = 17000, verbose = TRUE, FUN = Chisquare, ...)
}
\arguments{
  \item{x}{A vector of genotype counts (AA,AB,BB)}
  \item{nperm}{The number of permutations}
  \item{verbose}{\code{verbose} = TRUE will print results, \code{verbose} = FALSE is silent.}    
  \item{FUN}{An function call for calculating the test statistic for HWE
    (see examples below)}
  \item{\dots}{Additional parameters for the function call argument \code{FUN}}
}
\details{
The set of alleles for the observed sample is permuted. Consequently, the test is
conditional on allele frequency.
}
\value{Returns the p-value of the test}
\references{
Ziegler, A. & K\"onig, I.R. (2006) A statistical approach to genetic
epidemiology. Wiley.
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{\code{\link{HWChisq}},\code{\link{HWExact}},\code{\link{HWLratio}}}
\examples{
x <- c(MM=298,MN=489,NN=213)
\dontrun{
HW.test <- HWPerm(x,nperm=10000,verbose=TRUE) # uses default chi-square statistic
HW.test <- HWPerm(x,nperm=10000,verbose=TRUE,function(z)
HWChisq(z)$chisq,cc=0.5) # uses chi-square statistic with continuity correction.
HW.test <- HWPerm(x,nperm=10000,verbose=TRUE,function(y) HWLratio(y)$G2)
# uses likelihood ratio statistic.
HWPerm(x,nperm=10000,verbose=TRUE,function(y) 1-HWExact(y)$pval) # uses
exact test p-value
}
}
\keyword{htest}

