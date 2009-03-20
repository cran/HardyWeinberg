\name{HWData}
\alias{HWData}
\title{Generating marker data in Hardy-Weinberg Equilibrium}
\description{
   HWData generates samples of genotypic counts under various
   schemes. It mainly uses sampling from the multinomial distribution given
   Hardy-Weinberg allele frequencies.
}
\usage{
HWData(n = 100, nm = 100, f = 0, p = NULL, pfixed = FALSE, exactequilibrium = FALSE, pdist = "runif", ...)
}
\arguments{
  \item{n}{the sample size.}  
  \item{nm}{the number of markers (or samples).}  
  \item{f}{the inbreeding coefficient}
  \item{p}{the allele frequency}
  \item{pfixed}{if \code{TRUE} Haldane's distribution is used for sampling,
    if \code{FALSE} a multinomial distribution is used}
  \item{exactequilibrium}{generates data in exact HWE if set to \code{TRUE}}
  \item{pdist}{take a random allele frequency from a uniform or beta
    distribution of \code{pfixed = FALSE} and \code{p} is not given.}  
  \item{\dots}{specific parameters for the uniform or beta}  
}
\value{
  \item{Xt}{the genotypic counts.}
  \item{Xc}{the genotypic compositions.}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{HWTernaryPlot}
\examples{
n <- 100
nm <- 100
out <- HWData(n,nm)
}
\keyword{datagen}
