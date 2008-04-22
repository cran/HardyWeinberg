\name{HWData}
\alias{HWData}
\title{Generating marker data in Hardy-Weinberg Equilibrium}
\description{
   HWData takes samples from the multinomial distribution given
   Hardy-Weinberg allele frequencies.
}
\usage{
HWData(n = 100, nm = 100, pfixed = NULL, exactequilibrium = FALSE, pdist
= "runif",...)
}
\arguments{
  \item{n}{the sample size.}
  \item{nm}{the number of markers (or samples).}
  \item{pfixed}{take a fixed allele frequency with value \code{pfixed}.}
  \item{exactequilibrium}{generate data in exact HWE or use the multinomial distribution}
  \item{pdist}{take a random allele frequency from a uniform or beta
    distribution.}
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

