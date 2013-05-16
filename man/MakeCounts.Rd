\name{MakeCounts}
\alias{MakeCounts}

\title{Create genotype counts from bi-allelic marker data}
\description{
\code{MakeCounts} creates a matrix of genotype counts, with one row for
each bi-allelic marker, containing 4 columns with the counts AA, AB, BB
and NA (missings) respectively
}
\usage{
MakeCounts(X, alleles, pos1 = 1, pos2 = 3)
}
\arguments{
  \item{X}{A matrix or dataframe with bi-allelic genotyping information,
    markers in columns, individuals in rows}
  \item{alleles}{a vector of alleles for each marker (e.g. c("A/T","A/G",...))}
  \item{pos1}{position of the first allele in the allele string (1 by default)}
  \item{pos2}{position of the second allele in the allele string (3 by default)}
}
\details{
  \code{MakeCounts} is thought for bi-allelic marker data only. Missings
  are should be coded by NA. It produces the right input for \code{HWTernaryPlot}.
}
\value{A matrix of 4 columns}
\author{ Jan Graffelman \email{jan.graffelman@upc.edu} }
\seealso{\code{\link{HWTernaryPlot}}}
\examples{
SNP1 <- c("GG","GG","GG","GG","GG","GG","GG","GG","GG")
SNP2 <- c("CG","GG","CC","GG","GG","CG","CG","CG","CG")
SNP3 <- c("AA","AA","AA","AG","AA","AG","AA","AA","AA")
SNP4 <- c("GG","GG","GG","GG","GG","GG","GG","GG","GG")
SNP5 <- c("CC","CC","CC","CC","CC","CC","CT","CT","CT")
X <- cbind(SNP1,SNP2,SNP3,SNP4,SNP5)
Y <- MakeCounts(X,c("A/G","C/G","A/G","A/G","C/T"))
print(Y)
}
\keyword{misc}

