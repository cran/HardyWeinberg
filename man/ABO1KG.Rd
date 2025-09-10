\name{ABO1KG}
\alias{ABO1KG}
\docType{data}
\title{
Predicted ABO genotypes for the Expanded 1000 Genomes Project cohort
}
\description{
   The dataframe contains demographic variables, genotype data of five genetic polymorphisms and
   their predicted ABO blood group for 3202 individuals of the expanded 1000 Genomes Project cohort.
}
\usage{data("ABO1KG")}
\format{
  A data frame with 3202 observations on the following 12 variables.
  \describe{
    \item{\code{SAMPLE_NAME}}{The identifier for the individual of the 1KG project}
    \item{\code{POPULATION}}{Abbreviation of the population according to the 1KG project}
    \item{\code{fatherID}}{The identifier of the father of the individual (0 if missing)}
    \item{\code{motherID}}{The identifier of the mother of the individual (0 if missing)}
    \item{\code{sex}}{Sex (1=male; 2=female)}
    \item{\code{CONTINENT}}{Abbreviation of the continental region according to the 1KG project}
    \item{\code{rs8176719}}{a character vector}
    \item{\code{rs7853989}}{a character vector}
    \item{\code{rs8176743}}{a character vector}
    \item{\code{rs8176746}}{a character vector}
    \item{\code{rs8176747}}{a character vector}
    \item{\code{ABO}}{Predicted blood group}
  }
}
\details{
   The five polymorphisms have been used to create haplotypes which were used to assign the blood groups. Individuals with
   haplotypes corresponding to unknown bloodgroups were assigned a missing values (NA).
}
\source{
   https://www.internationalgenome.org/
}
\references{
   Byrska-Bishop, M. et al. (2022) High-coverage whole-genome sequencing of the expanded 1000 Genomes Project cohort including 602 trios. Cell 185(18) pp. 3426-3440.e19. doi 10.1016/j.cell.2022.08.004.

   Graffelman, J. (2025) A logratio approach to the analysis of autosomal genotype frequencies across multiple samples. Preprint available at www.bioRxiv.org. doi: 10.1101/2025.01.17.633675
}
\examples{
  data(ABO1KG)
  str(ABO1KG)
}
\keyword{datasets}
