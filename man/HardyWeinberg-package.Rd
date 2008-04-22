\name{HardyWeinberg-package}
\alias{HardyWeinberg-package}
\alias{HardyWeinberg}
\docType{package}
\title{Graphical tests for Hardy-Weinberg equilibrium}

\description{
The HardyWeinberg-package provides graphical tests for Hardy-Weinberg equilibrium (HWE) based 
on the ternary plot (de Finetti diagram). The package constructs ternary plots for genotypic 
compositions for bi-allelic marker data. The acceptance region for several statistical tests
of HWE (Chisquare, Chisquare with continuity correction, Fisher exact) can be depicted inside
the ternary plot with the routines of the package. Large numbers of bi-allelic markers 
(e.g. SNPs) can be represented in a single ternary diagram and the statistical (non)significance
of a test for HWE can be inferred from the position of the marker in the plot.
}
\details{
\tabular{ll}{
Package: \tab HardyWeinberg\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2008-04-22\cr
License: \tab GPL Version 2 or later.\cr
}

The most important function of the package is \code{HWTernaryPlot} that can
be used to create ternary plots with acceptance regions for HWE. Other routines
implement statistical tests for HWE such as \code{HWChisq} and \code{HWLratio}.

}
\author{
Jan Graffelman

Maintainer: Jan Graffelman <jan.graffelman@upc.edu>

}
\references{
Weir, B.S. (1996) \emph{Genetic Data Analysis II}. Sinauer Associates, Massachusetts.

Graffelman, J. and Morales, J. (2008) Graphical tests for Hardy-Weinberg equilibrium 
based on the ternary plot. \emph{Human Heredity} 65(2):77-84.
}
\keyword{ package }
\examples{

library(HardyWeinberg)

# make random compositions that are in HWE

set.seed(123)

m <- 100 # number of markers
n <- 100 # sample size

Xc <- HWData(n,m)$Xc
out <- HWTernaryPlot(Xc,100,region=1,vertex.cex=2,signifcolour=TRUE)

}
