\name{TSIXTriAllelics}
\docType{data}
\alias{TSIXTriAllelics}
\title{
  Tri-allelic polymorphisms on the X chromosome of the TSI population
}
\description{
  This dataframe contains genotype counts for six three-allelic polymorphisms (A,B,C) on chromosome X of a sample of individuals from the TSI population (Tuscany, Italy) of the 1,000 genomes project.
}
\usage{data(TSIXTriAllelics)}
\format{
  A data frame with 6 observations on the following 10 variables.
  \describe{
    \item{\code{ID}}{Identifier of the polymorphism}
    \item{\code{A}}{Male A genotype count}
    \item{\code{B}}{Male B genotype count}
    \item{\code{C}}{Male C genotype count}
    \item{\code{AA}}{Female AA genotype count}
    \item{\code{AB}}{Female AB genotype count}
    \item{\code{AC}}{Female AC genotype count}
    \item{\code{BB}}{Female BB genotype count}
    \item{\code{BC}}{Female BC genotype count}
    \item{\code{CC}}{Female CC genotype count}
  }
}
\source{
  Data taken from the 1,000 genomes project at www.internationalgenome.org
}
\references{
  Graffelman, J. and Ortoleva, L. (2020) A network algorithm for the X chromosomal exact test for Hardy-Weinberg equilibrium with multiple alleles. Under review.
}
\examples{
data(TSIXTriAllelics)
}
\keyword{datasets}
