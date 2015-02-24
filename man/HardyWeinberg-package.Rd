\name{HardyWeinberg-package}
\alias{HardyWeinberg-package}
\alias{HardyWeinberg}
\docType{package}
\title{Graphical tests for Hardy-Weinberg equilibrium}

\description{
 The package HardyWeinberg offers tools for exploring
 diallelic genetic marker data. It offers all classical tests (chi-square,
 exact, likelihood-ratio and permutation tests) for
 Hardy-Weinberg equilibrium, functions for power computation and
 for the simulation of marker data under equilibrium and disequilibrium.
 Functions for testing equilibrium in the presence of missing data by
 using multiple imputation are provided. The package also supplies
 various graphical tools such as ternary plots with acceptance regions,
 log-ratio plots and Q-Q plots for exploring the equilibrium status of
 a large set of diallelic markers. 
}
\details{
\tabular{ll}{
Package: \tab HardyWeinberg\cr
Type: \tab Package\cr
Version: \tab 1.5.5\cr
Date: \tab 2015-02-24\cr
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
based on the ternary plot. \emph{Human Heredity} 65(2):77-84. \url{http://dx.doi.org/10.1159/000108939}.

Graffelman, J. (2015) Exploring Diallelic Genetic Markers: The {HardyWeinberg} Package.
\emph{Journal of Statistical Software} 64(3): 1-22. \url{http://www.jstatsoft.org/v64/i03/}.
}

\keyword{ package }
\examples{

library(HardyWeinberg)

# make random compositions that are in HWE

set.seed(123)

m <- 100 # number of markers
n <- 100 # sample size

out <- HWData(n,m)
Xc <- out$Xc
out <- HWTernaryPlot(Xc,100,region=1,vertex.cex=2,signifcolour=TRUE)

}
