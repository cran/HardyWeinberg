## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.align = "center")
knitr::opts_chunk$set(warnings = FALSE)
knitr::opts_chunk$set(message = FALSE)
knitr::opts_chunk$set(fig.width = 6, fig.height = 6) 

## ----preinstall---------------------------------------------------------------
#install.packages("HardyWeinberg")
library(HardyWeinberg)

## -----------------------------------------------------------------------------
x <- c(MM = 298, MN = 489, NN = 213)

## -----------------------------------------------------------------------------
maf(x)

## -----------------------------------------------------------------------------
maf(x,2)
maf(x,3)

## -----------------------------------------------------------------------------
HWTernaryPlot(x,region=0,hwcurve=FALSE,grid=TRUE,markercol="blue")

## -----------------------------------------------------------------------------
HW.test <- HWChisq(x, verbose = TRUE)

## -----------------------------------------------------------------------------
HW.test <- HWChisq(x, cc = 0, verbose = TRUE)

## -----------------------------------------------------------------------------
HW.lrtest <- HWLratio(x, verbose = TRUE)

## -----------------------------------------------------------------------------
HW.exacttest <- HWExact(x, verbose = TRUE, pvaluetype = "midp")

## -----------------------------------------------------------------------------
set.seed(123)
#HW.permutationtest <- HWPerm(x, verbose = TRUE)
#Permutation test for Hardy-Weinberg equilibrium
#Observed statistic: 0.2214896   17000 permutations. p-value: 0.6508235 

## ----alltests-----------------------------------------------------------------
#HW.results <- HWAlltests(x, verbose = TRUE, include.permutation.test = TRUE)
#                                            Statistic   p-value
#Chi-square test:                            0.2214896 0.6379073
#Chi-square test with continuity correction: 0.1789563 0.6722717
#Likelihood-ratio test:                      0.2214663 0.6379250
#Exact test with selome p-value:                    NA 0.6556635
#Exact test with dost p-value:                      NA 0.6723356
#Exact test with mid p-value:                       NA 0.6330965
#Permutation test:                           0.2214896 0.6508235

## -----------------------------------------------------------------------------
post.dens <- HWLindley(seq(-1,1,by=0.01),x)
plot(seq(-1,1,by=0.01),post.dens,type="l",xlab=expression(alpha),
     ylab=expression(pi(alpha)))
segments(0,0,0,HWLindley(0,x),lty="dotted",col="red")

## -----------------------------------------------------------------------------
HWLindley.cri(x=x)

## -----------------------------------------------------------------------------
SNP1 <- c(A=399,B=205,AA=230,AB=314,BB=107) 

## -----------------------------------------------------------------------------
HWChisq(SNP1,cc=0,x.linked=TRUE,verbose=TRUE)

## -----------------------------------------------------------------------------
HWChisq(SNP1[3:5],cc=0)

## -----------------------------------------------------------------------------
HWLratio(SNP1,x.linked=TRUE)

## -----------------------------------------------------------------------------
HWExact(SNP1,x.linked=TRUE)

## -----------------------------------------------------------------------------
HWExact(SNP1,x.linked=TRUE,pvaluetype="midp")

## -----------------------------------------------------------------------------
HWExact(SNP1[3:5],pvaluetype="midp")

## ----permutationxlinked-------------------------------------------------------
#HWPerm(SNP1,x.linked=TRUE)
#Permutation test for Hardy-Weinberg equilibrium of an X-linked marker
#Observed statistic: 7.624175   17000 permutations. p-value: 0.02211765 

## ----alltestsxlinked----------------------------------------------------------
#HWAlltests(SNP1,x.linked=TRUE,include.permutation.test=TRUE)
#                                            Statistic    p-value
#Chi-square test:                             7.624175 0.02210200
#Chi-square test with continuity correction:  7.242011 0.02675576
#Likelihood-ratio test:                       7.693436 0.02134969
#Exact test with selome p-value:                    NA 0.02085798
#Exact test with dost p-value:                      NA         NA
#Exact test with mid p-value:                       NA 0.02082957
#Permutation test:                            7.624175 0.02211765

## -----------------------------------------------------------------------------
AFtest(SNP1)

## ----bayesianx----------------------------------------------------------------
HWPosterior(males=SNP1[1:2],females=SNP1[3:5],x.linked=TRUE)

## -----------------------------------------------------------------------------
data(Markers)
Markers[1:12,]

## -----------------------------------------------------------------------------
Xt <- table(Markers[,1])
Xv <- as.vector(Xt)
names(Xv) <- names(Xt)
Xv
HW.test <- HWChisq(Xv,cc=0,verbose=TRUE)

## -----------------------------------------------------------------------------
set.seed(123)
Results <- HWMissing(Markers[,1], m = 50, method = "sample", verbose=TRUE)

## -----------------------------------------------------------------------------
set.seed(123)
Results <- HWMissing(Markers[, 1:5], m = 50, verbose = TRUE)

## -----------------------------------------------------------------------------
set.seed(123)
Results <- HWMissing(Markers[, 1:5], m = 50, statistic = "exact", verbose = TRUE)

## -----------------------------------------------------------------------------
x <- c(MM = 298, MN = 489, NN = 213)
n <- sum(x)
nM <- mac(x) 
pw4 <- HWPower(n, nM, alpha = 0.05, test = "exact", theta = 4, 
               pvaluetype = "selome")
print(pw4)
pw8 <- HWPower(n, nM, alpha = 0.05, test = "exact", theta = 8, 
               pvaluetype = "selome")
print(pw8)

## ----glyoxalasa---------------------------------------------------------------
data("Glyoxalase")
Glyoxalase <- as.matrix(Glyoxalase)
HWStrata(Glyoxalase)

## -----------------------------------------------------------------------------
pvalues <- HWExactStats(Glyoxalase)
pvalues

## -----------------------------------------------------------------------------
g1 <- c(AA=0.034, AB=0.330, BB=0.636)
g2 <- c(AA=0.349, AB=0.493, BB=0.158)
x  <- c(AA=0.270, AB=0.453, BB=0.277)

## -----------------------------------------------------------------------------
G <- cbind(g1,g2)
contributions <- HWEM(x,G=G)
contributions

## -----------------------------------------------------------------------------
p <- c(af(g1),af(g2))
contributions <- HWEM(x,p=p)
contributions

## -----------------------------------------------------------------------------
data(JPTsnps)
JPTsnps[1,]
Results <- HWPosterior(males=JPTsnps[1,1:3],
                       females=JPTsnps[1,4:6],
                       x.linked=FALSE,precision=0.05)

## -----------------------------------------------------------------------------
data(JPTsnps)
AICs <- HWAIC(JPTsnps[1,1:3],JPTsnps[1,4:6])
AICs

## ----datasets-----------------------------------------------------------------
data("CEUchr22")
CEUchr22[1:5,1:5]

## -----------------------------------------------------------------------------
Z <- MakeCounts(CEUchr22)
Z <- Z[,1:3]
head(Z)

## -----------------------------------------------------------------------------
HWTernaryPlot(Z[,1:3],patternramp=TRUE,region=0)

## -----------------------------------------------------------------------------
pminor <- maf(Z)
HWTernaryPlot(Z[pminor>0.05,],patternramp=TRUE,region=0)

## -----------------------------------------------------------------------------
hist(pminor[pminor > 0],freq=FALSE,xlab="MAF",main="CEU MAF CHR 22")

## -----------------------------------------------------------------------------
cminor <- maf(Z[pminor > 0,],option=3)
barplot(table(cminor[,1]),cex.names = 0.75)

## ----manytests----------------------------------------------------------------
Zpoly <- Z[!is.mono(Z),]
npoly <- nrow(Zpoly)
chisq.pvalues <- HWChisqStats(Zpoly,pvalues=TRUE)
exact.pvalues <- HWExactStats(Zpoly,midp=TRUE)
bonferronithreshold <- 0.05/npoly
sum(chisq.pvalues < bonferronithreshold)
sum(exact.pvalues < bonferronithreshold)

## ----ternaryCEU---------------------------------------------------------------
Zu <- UniqueGenotypeCounts(Z)
Zu <- Zu[,1:3]
HWTernaryPlot(Zu,region=1,pch=1)

## -----------------------------------------------------------------------------
HWTernaryPlot(Zu[,1:3],region=7,pch=1)

## -----------------------------------------------------------------------------
data("CEUchr22")

Z <- MakeCounts(CEUchr22)
Z <- Z[,1:3]
Z <- Z[!is.mono(Z),]

alfreq  <- af(Z)
alcount <- maf(Z,option=3)[,1]

chisq.pvals <- HWChisqStats(Z,pvalues=TRUE)

set.seed(123)
Z.sim.chi <- HWData(nm=nrow(Z),n=99,p=alfreq)
Z.sim.chi <- Z.sim.chi[!is.mono(Z.sim.chi),]
chisq.pvals.sim <- HWChisqStats(Z.sim.chi,pvalues=TRUE)

set.seed(123)
Z.sim.exa <- HWData(nm=nrow(Z),n=99,nA=alcount,conditional = TRUE)
Z.sim.exa <- Z.sim.exa[!is.mono(Z.sim.exa),]

opar <- par(mfrow=c(2,2),mar=c(3,3,2,0)+0.5,mgp=c(2,1,0))
par(mfg=c(1,1))
 qqunif(chisq.pvals,logplot = TRUE,
       plotline = 1,main="A: observed QQ uniform")
par(mfg=c(1,2))
 qqunif(chisq.pvals.sim,logplot = TRUE,xylim=15,
        main="B: simulated QQ uniform")
par(mfg=c(2,1))
 set.seed(123)
 HWQqplot(Z,nsim=5,logplot=TRUE,main="C: observed QQ HWE null")
par(mfg=c(2,2))
 set.seed(123)
 HWQqplot(Z.sim.exa,nsim=5,logplot=TRUE,main="D: simulated QQ HWE null")
par(opar)

## ----simulate-----------------------------------------------------------------
set.seed(123)
n <- 100
m <- 100
X1 <- HWData(m, n, p = rep(0.5, m))
X2 <- HWData(m, n)
X3 <- HWData(m, n, p = rep(0.5, m), f = rep(0.5, m))
X4 <- HWData(m, n, f = rep(0.5, m))
X5 <- HWData(m, n, p = rep(c(0.2, 0.4, 0.6, 0.8), 25), conditional = TRUE)
X6 <- HWData(m, n, exactequilibrium = TRUE)
opar <- par(mfrow = c(3, 2),mar = c(1, 0, 3, 0) + 0.1)
par(mfg = c(1, 1))
HWTernaryPlot(X1, main = "(a)")
par(mfg = c(1, 2))
HWTernaryPlot(X2, main = "(b)")
par(mfg = c(2, 1))
HWTernaryPlot(X3, main = "(c)")
par(mfg = c(2, 2))
HWTernaryPlot(X4, main = "(d)")
par(mfg = c(3, 1))
HWTernaryPlot(X5, main = "(e)")
par(mfg = c(3, 2))
HWTernaryPlot(X6, main = "(f)")
par(opar)

## -----------------------------------------------------------------------------
X <- HWData(nm=100,shape1=1,shape2=10)

## -----------------------------------------------------------------------------
x <- c(fA=182,fB=60,nAB=17,nOO=176)
al.fre <- HWABO(x)

## -----------------------------------------------------------------------------
x <- c(AA=20,AB=31,AC=26,BB=15,BC=12,CC=0)
#results <- HWTriExact(x)
#Tri-allelic Exact test for HWE (autosomal).
#Allele counts: A = 38 B = 73 C = 97 
#sum probabilities all outcomes 1 
#probability of the sample 0.0001122091 
#p-value =  0.03370688 

## -----------------------------------------------------------------------------
x <- c(AA=20,AB=31,AC=26,BB=15,BC=12,CC=0)
x <- toTriangular(x)
x
m <- c(A=0,B=0,C=0)
results <- HWNetwork(ma=m,fe=x)

## -----------------------------------------------------------------------------
set.seed(123)
x <- c(AA=20,AB=31,AC=26,BB=15,BC=12,CC=0)
x <- toTriangular(x)
#results <- HWPerm.mult(x)
#Permutation test for Hardy-Weinberg equilibrium (autosomal).
#3 alleles detected.
#Observed statistic: 0.0001122091   17000 permutations. p-value: 0.03405882 

## -----------------------------------------------------------------------------
males   <- c(A=1,B=21,C=34) 
females <- c(AA=0,AB=1,AC=0,BB=8,BC=24,CC=15)
results <- HWTriExact(females,males)

## -----------------------------------------------------------------------------
males   <- c(A=1,B=21,C=34) 
females <- toTriangular(c(AA=0,AB=1,AC=0,BB=8,BC=24,CC=15))
results <- HWNetwork(ma=males,fe=females)

## -----------------------------------------------------------------------------
data("NistSTRs")
NistSTRs[1:5,1:5]
n <- nrow(NistSTRs)
p <- ncol(NistSTRs)/2
n
p

A1 <- NistSTRs[,1]
A2 <- NistSTRs[,2]
GenotypeCounts <- AllelesToTriangular(A1,A2)
print(GenotypeCounts)

set.seed(123)
#out <- HWPerm.mult(GenotypeCounts)
#Permutation test for Hardy-Weinberg equilibrium (autosomal).
#7 alleles detected.
#Observed statistic: 2.290724e-11   17000 permutations. p-value: 0.8644706 

## -----------------------------------------------------------------------------
#Results <- HWStr(NistSTRs,test="permutation")
#29 STRs detected.
#> Results
#          STR   N Nt MinorAF MajorAF     Ho     He     Hp   pval
#1    CSF1PO-1 361  7  0.0055  0.3601 0.7341 0.7194 1.4016 0.8614
#2  D10S1248-1 361  9  0.0014  0.3075 0.7645 0.7586 1.5552 0.6488
#3   D12S391-1 361 16  0.0014  0.1717 0.8975 0.8909 2.3686 0.8955
#4   D13S317-1 361  8  0.0014  0.3255 0.7618 0.7837 1.7102 0.1049
#5   D16S539-1 361  7  0.0180  0.3144 0.7645 0.7600 1.5933 0.1641
#6    D18S51-1 361 15  0.0014  0.1704 0.8560 0.8758 2.2139 0.1926
#7   D19S433-1 361 15  0.0014  0.3615 0.7673 0.7694 1.7624 0.9667
#8   D1S1656-1 361 15  0.0028  0.1496 0.9252 0.8992 2.4073 0.8699
#9    D21S11-1 361 16  0.0014  0.2825 0.8227 0.8288 1.9946 0.6156
#10 D22S1045-1 361  8  0.0055  0.3823 0.7535 0.7220 1.4823 0.9785
#11  D2S1338-1 361 12  0.0014  0.1856 0.8698 0.8814 2.2464 0.5066
#12   D2S441-1 361 11  0.0014  0.3435 0.7867 0.7693 1.6734 0.7125
#13  D3S1358-1 361  9  0.0014  0.2729 0.7562 0.7900 1.6437 0.3789
#14   D5S818-1 361  9  0.0014  0.3878 0.7064 0.6977 1.3939 0.7615
#15  D6S1043-1 361 14  0.0014  0.2964 0.7978 0.8232 2.0059 0.0220
#16   D7S820-1 361  9  0.0014  0.2562 0.8310 0.8161 1.7925 0.5045
#17  D8S1179-1 361 10  0.0014  0.3296 0.7839 0.8072 1.8386 0.8606
#18   F13A01-1 361 12  0.0014  0.3490 0.7590 0.7353 1.5471 0.8793
#19     F13B-1 361  6  0.0055  0.3892 0.7008 0.7175 1.3790 0.8276
#20   FESFPS-1 361  7  0.0014  0.4114 0.6731 0.6930 1.3085 0.8506
#21      FGA-1 361 14  0.0014  0.2050 0.8670 0.8594 2.1163 0.0528
#22      LPL-1 361  8  0.0014  0.4224 0.7175 0.6947 1.3386 0.9721
#23  Penta_C-1 361 10  0.0014  0.3947 0.7701 0.7523 1.6035 0.0248
#24  Penta_D-1 361 13  0.0014  0.2327 0.8587 0.8247 1.8925 0.1464
#25  Penta_E-1 361 19  0.0014  0.1994 0.8920 0.8910 2.4391 0.4811
#26     SE33-1 361 39  0.0014  0.0942 0.9501 0.9483 3.1472 0.2506
#27     TH01-1 361  8  0.0014  0.3449 0.7424 0.7646 1.5616 0.2379
#28     TPOX-1 361  8  0.0014  0.5249 0.6537 0.6404 1.2572 0.5545
#29      vWA-1 361 10  0.0014  0.2839 0.8061 0.8076 1.7577 0.4501

