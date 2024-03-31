### R code from vignette source 'HardyWeinberg.Rnw'

###################################################
### code chunk number 1: HardyWeinberg.Rnw:604-605
###################################################
options(prompt = "R> ", continue = "+ ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: HardyWeinberg.Rnw:608-610 (eval = FALSE)
###################################################
## install.packages("HardyWeinberg")
## library("HardyWeinberg")


###################################################
### code chunk number 3: HardyWeinberg.Rnw:618-619 (eval = FALSE)
###################################################
## vignette("HardyWeinberg")


###################################################
### code chunk number 4: HardyWeinberg.Rnw:629-631
###################################################
library("HardyWeinberg")
x <- c(MM = 298, MN = 489, NN = 213)


###################################################
### code chunk number 5: HardyWeinberg.Rnw:638-639
###################################################
HW.test <- HWChisq(x, verbose = TRUE)


###################################################
### code chunk number 6: HardyWeinberg.Rnw:645-646
###################################################
HW.test <- HWChisq(x, cc = 0, verbose = TRUE)


###################################################
### code chunk number 7: HardyWeinberg.Rnw:656-657
###################################################
HW.lrtest <- HWLratio(x, verbose = TRUE)


###################################################
### code chunk number 8: HardyWeinberg.Rnw:668-669
###################################################
HW.exacttest <- HWExact(x, verbose = TRUE)


###################################################
### code chunk number 9: HardyWeinberg.Rnw:679-681
###################################################
set.seed(123)
HW.permutationtest <- HWPerm(x, verbose = TRUE)


###################################################
### code chunk number 10: HardyWeinberg.Rnw:688-690
###################################################
x <- c(MN = 489, NN = 213, MM = 298)
HW.test <- HWChisq(x, verbose = TRUE)


###################################################
### code chunk number 11: HardyWeinberg.Rnw:698-699
###################################################
HW.results <- HWAlltests(x, verbose = TRUE, include.permutation.test = TRUE)


###################################################
### code chunk number 12: HardyWeinberg.Rnw:708-710
###################################################
data(Markers)
Markers[1:12,]


###################################################
### code chunk number 13: HardyWeinberg.Rnw:716-720
###################################################
Xt <- table(Markers[,1])
Xv <- as.vector(Xt)
names(Xv) <- names(Xt)
HW.test <- HWChisq(Xv,cc=0,verbose=TRUE)


###################################################
### code chunk number 14: HardyWeinberg.Rnw:727-729
###################################################
set.seed(123)
Results <- HWMissing(Markers[,1], m = 50, method = "sample", verbose=TRUE)


###################################################
### code chunk number 15: HardyWeinberg.Rnw:736-738
###################################################
set.seed(123)
Results <- HWMissing(Markers[, 1:5], m = 50, verbose = TRUE)


###################################################
### code chunk number 16: HardyWeinberg.Rnw:744-746
###################################################
set.seed(123)
Results <- HWMissing(Markers[, 1:5], m = 50, statistic = "exact", verbose = TRUE)


###################################################
### code chunk number 17: HardyWeinberg.Rnw:755-757
###################################################
data(JPTsnps)
Results <- HWPosterior(males=JPTsnps[1,1:3],females=JPTsnps[1,4:6],x.linked=FALSE,precision=0.05)


###################################################
### code chunk number 18: HardyWeinberg.Rnw:764-767
###################################################
data(JPTsnps)
AICs <- HWAIC(JPTsnps[1,1:3],JPTsnps[1,4:6])
AICs


###################################################
### code chunk number 19: HardyWeinberg.Rnw:777-780
###################################################
g1 <- c(0.034, 0.330, 0.636)
g2 <- c(0.349, 0.493, 0.158)
x  <- c(0.270, 0.453,0.277)


###################################################
### code chunk number 20: HardyWeinberg.Rnw:785-788
###################################################
G <- cbind(g1,g2)
contributions <- HWEM(x,G=G)
contributions


###################################################
### code chunk number 21: HardyWeinberg.Rnw:793-796
###################################################
p <- c(af(g1),af(g2))
contributions <- HWEM(x,p=p)
contributions


###################################################
### code chunk number 22: HardyWeinberg.Rnw:805-807
###################################################
SNP1 <- c(A=399,B=205,AA=230,AB=314,BB=107) 
HWChisq(SNP1,cc=0,x.linked=TRUE,verbose=TRUE)


###################################################
### code chunk number 23: HardyWeinberg.Rnw:812-813
###################################################
HWChisq(SNP1[3:5],cc=0)


###################################################
### code chunk number 24: HardyWeinberg.Rnw:821-822
###################################################
HWExact(SNP1,x.linked=TRUE)


###################################################
### code chunk number 25: HardyWeinberg.Rnw:827-828
###################################################
HWExact(SNP1,x.linked=TRUE,pvaluetype="midp")


###################################################
### code chunk number 26: HardyWeinberg.Rnw:834-835
###################################################
HWExact(SNP1[3:5])


###################################################
### code chunk number 27: HardyWeinberg.Rnw:840-841
###################################################
HWPerm(SNP1,x.linked=TRUE)


###################################################
### code chunk number 28: HardyWeinberg.Rnw:846-847
###################################################
HWLratio(SNP1,x.linked=TRUE)


###################################################
### code chunk number 29: HardyWeinberg.Rnw:852-853
###################################################
HWAlltests(SNP1,x.linked=TRUE,include.permutation.test=TRUE)


###################################################
### code chunk number 30: HardyWeinberg.Rnw:858-859
###################################################
AFtest(SNP1)


###################################################
### code chunk number 31: HardyWeinberg.Rnw:869-870
###################################################
HWPosterior(males=SNP1[1:2],females=SNP1[3:5],x.linked=TRUE)


###################################################
### code chunk number 32: HardyWeinberg.Rnw:889-898
###################################################
x <- c(MM = 298, MN = 489, NN = 213)
n <- sum(x)
nM <- mac(x) 
pw4 <- HWPower(n, nM, alpha = 0.05, test = "exact", theta = 4, 
               pvaluetype = "selome")
print(pw4)
pw8 <- HWPower(n, nM, alpha = 0.05, test = "exact", theta = 8, 
               pvaluetype = "selome")
print(pw8)


###################################################
### code chunk number 33: HardyWeinberg.Rnw:917-940
###################################################
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


###################################################
### code chunk number 34: HardyWeinberg.Rnw:946-947
###################################################
X <- HWData(nm=100,shape1=1,shape2=10)


###################################################
### code chunk number 35: HardyWeinberg.Rnw:955-958 (eval = FALSE)
###################################################
## data("HapMapCHBChr1", package = "HardyWeinberg")
## HWTernaryPlot(HapMapCHBChr1, region = 1)
## HWTernaryPlot(HapMapCHBChr1, region = 7)


###################################################
### code chunk number 36: HardyWeinberg.Rnw:974-981 (eval = FALSE)
###################################################
## set.seed(123)
## data("HapMapCHBChr1", package = "HardyWeinberg")
## HWQqplot(HapMapCHBChr1)
## dev.off()
## set.seed(123)
## SimulatedData <- HWData(nm = 225, n = 84, p = af(HapMapCHBChr1))
## HWQqplot(SimulatedData)


###################################################
### code chunk number 37: HardyWeinberg.Rnw:1004-1006
###################################################
x <- c(fA=182,fB=60,nAB=17,nOO=176)
al.fre <- HWABO(x)


###################################################
### code chunk number 38: HardyWeinberg.Rnw:1015-1017
###################################################
x <- c(AA=20,AB=31,AC=26,BB=15,BC=12,CC=0)
results <- HWTriExact(x)


###################################################
### code chunk number 39: HardyWeinberg.Rnw:1022-1026
###################################################
x <- c(AA=20,AB=31,AC=26,BB=15,BC=12,CC=0)
x <- toTriangular(x)
m <- c(A=0,B=0,C=0)
#results <- HWNetwork(ma=m,fe=x)


###################################################
### code chunk number 40: HardyWeinberg.Rnw:1031-1034
###################################################
males   <- c(A=1,B=21,C=34) 
females <- c(AA=0,AB=1,AC=0,BB=8,BC=24,CC=15)
results <- HWTriExact(females,males)


###################################################
### code chunk number 41: HardyWeinberg.Rnw:1039-1042
###################################################
males   <- c(A=1,B=21,C=34) 
females <- toTriangular(c(AA=0,AB=1,AC=0,BB=8,BC=24,CC=15))
#results <- HWNetwork(ma=males,fe=females)


###################################################
### code chunk number 42: HardyWeinberg.Rnw:1047-1051
###################################################
set.seed(123)
x <- c(AA=20,AB=31,AC=26,BB=15,BC=12,CC=0)
x <- toTriangular(x)
#results <- HWPerm.mult(x)


###################################################
### code chunk number 43: HardyWeinberg.Rnw:1058-1065
###################################################
set.seed(123)
data(NistSTRs)
A1 <- NistSTRs[,1]
A2 <- NistSTRs[,2]
GenotypeCounts <- AllelesToTriangular(A1,A2)
print(GenotypeCounts)
#out <- HWPerm.mult(GenotypeCounts)


###################################################
### code chunk number 44: HardyWeinberg.Rnw:1070-1071
###################################################
#Results <- HWStr(NistSTRs,test="permutation")


