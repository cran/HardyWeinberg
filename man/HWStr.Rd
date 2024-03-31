\name{HWStr}
\alias{HWStr}
\title{
Testing a set of microsatellites (STRs) for Hardy-Weinberg equilibrium 
}
\description{
Function \code{HWStr} is a wrapper function arount \code{HWPerm.mult} and \code{HWChisq} in order to test a set of STRs for Hardy-Weinberg equilibrium.
}
\usage{
HWStr(X, test = "permutation", verbose = FALSE, ...)
}
\arguments{
  \item{X}{A data matrix with STRs in columns, with there two alleles coded in successive columns.}
  \item{test}{\code{permutation} for a permutation test, or \code{"chisq"} for a chisquare test}
  \item{verbose}{be silent if set to \code{FALSE} or print info for each STR if \code{TRUE}.}
  \item{\dots}{possible extra parameters to be passed on to \code{HWPerm.mult} such as \code{nperm} or \code{eps}.}
}
\details{It is recommended to test by a permutation test. By default, 17.000 permutations are used. Exact testing is not implemented.
         Missing genotypes are excluded on a per STR basis.
  }
\value{A data frame with the 
\item{strnames}{STR name}
\item{N}{sample size}
\item{Nt}{number of alleles}
\item{MinorAF}{minor allele frequency}
\item{MajorAF}{major allele frequency}
\item{Ho}{observed heterozygosity}
\item{He}{expected heterozygosity (without bias correction)}
\item{Hp}{Shannon index of allele frequencies}
\item{pval}{p-value of the HWE test}
}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{
 \code{\link{HWPerm.mult}}, \code{\link{HWChisq}}
}
\examples{
data(NistSTRs)
\dontrun{
Results <- HWStr(NistSTRs,test="permutation")
}
}
\keyword{htest}
