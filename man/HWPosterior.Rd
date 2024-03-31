\name{HWPosterior}
\alias{HWPosterior}
\title{
Calculation of posterior probabilities and Bayes factors for
Hardy-Weinberg tests at X-chromosomal variants.
}
\description{
Function \code{HWPosterior} calculates posterior probabilities and Bayes
factors for tests for Hardy-Weinberg equilibrium of autosomal and X-chromosomal
variants.
}
\usage{
HWPosterior(males, females, verbose = TRUE, prior.af = c(0.5,0.5), prior.gf =
            c(0.333,0.333,0.333), x.linked = FALSE, precision = 0.05) 
}
\arguments{
  \item{males}{A vector of male genotype counts. For autosomal variants, this
    should be a three-element named vector like (AA,AB,BB); for X-chromosomal
    variants it should be a two-element vector giving the counts of the hemizygous
    genotypes like (A,B).}
  \item{females}{A vector of female genotype counts. This should be a three-element
    named vector like (AA,AB,BB)} 
  \item{verbose}{prints results if \code{verbose = TRUE}}
  \item{prior.af}{Beta prior parameters for male and female allele frequencies}
  \item{prior.gf}{Dirichlet prior parameters for female genotype
    frequencies}
  \item{x.linked}{logical indicating whether the variant is autosomal or
    X-chromosomal}
  \item{precision}{precision parameter for marginal likelihoods that
    require numeric integration}
}
\details{For X-chromosomal variants, four possible models are considered, and the posterior
  probabilities and Bayes factors for each model are calculated.

  For autosomal variants, ten possible scenarios are considered, and the
  posterior probabilities for all models are calculated.

  In general, default Dirichlet priors are used for genotype
  frequencies, and beta prior are used for allele frequencies.
}
\value{
  For X-chromosomal variants, a matrix with posterior probabilities and
  Bayes factors will be produced.
  For autosomal variants, a vector of posterior probabilities is produced.
}
\references{
  Puig, X., Ginebra, J. and Graffelman, J. (2017) A Bayesian test for Hardy-Weinberg
  equilibrium of bi-allelic X-chromosomal markers. Heredity 119(4):226--236.
  \doi{10.1038/hdy.2017.30}.

  Puig, X., Ginebra, J. and Graffelman, J. (2019) Bayesian model selection for the study
  of Hardy-Weinberg proportions and homogeneity of gender allele frequencies.
  Heredity  123(5), pp. 549-564. \doi{10.1038/s41437-019-0232-0}
}
\author{
  Xavi Puig \email{xavier.puig@upc.edu} and
  Jan Graffelman \email{jan.graffelman@upc.edu}
}
\seealso{ \code{\link{HWChisq}}, \code{\link{HWExact}}, \code{\link{HWExactStats}} }
\examples{
#
# An X-chromosomal example
#
males   <- c(A=43,B=13)
females <- c(AA=26,AB=19,BB=3)
out     <- HWPosterior(males,females,verbose=TRUE,x.linked=TRUE)
#
# An autosomal example
#
data(JPTsnps)
males     <- JPTsnps[1,1:3]
females   <- JPTsnps[1,4:6]
post.prob <- HWPosterior(males,females,x.linked=FALSE)
}
\keyword{htest}

