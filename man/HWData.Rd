\name{HWData}
\alias{HWData}
\title{Generate genetic marker data in or out of Hardy-Weinberg Equilibrium}
\description{
   HWData generates samples of genotypic counts under various
   schemes. It mainly uses sampling from the multinomial distribution given
   Hardy-Weinberg allele frequencies.
}
\usage{
HWData(nm = 100, n = rep(100, nm), f = rep(0, nm), p = runif(nm),
       pfixed = FALSE, exactequilibrium = FALSE, pdist = "runif", ...)
}
\arguments{
  \item{nm}{the number of markers (or samples).}
  \item{n}{the sample size.}    
  \item{f}{the inbreeding coefficient}
  \item{p}{a vector of allele frequencies}
  \item{pfixed}{if \code{TRUE} Haldane's distribution is used for sampling,
    if \code{FALSE} a multinomial distribution is used}
  \item{exactequilibrium}{generates data in exact HWE if set to
    \code{TRUE}}
  \item{pdist}{take a random allele frequency from a uniform or beta
    distribution of \code{pfixed = FALSE} and \code{p} is not given.}
  \item{\dots}{specific parameters for the uniform or beta}  
}
\value{
  \item{Xt}{the genotypic counts.}
  \item{Xc}{the genotypic compositions.}
}
\author{ Jan Graffelman (jan.graffelman@upc.edu) }
\seealso{\code{\link{HWTernaryPlot}}}
\examples{
nm <- 100
n <- 100
out <- HWData(nm,n)
}
\keyword{datagen}
