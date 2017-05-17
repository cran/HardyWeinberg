### R code from vignette source 'HardyWeinberg.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: HardyWeinberg.Rnw:739-740
###################################################
options(prompt = "R> ", continue = "+ ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: HardyWeinberg.Rnw:743-745 (eval = FALSE)
###################################################
## install.packages("HardyWeinberg")
## library("HardyWeinberg")


###################################################
### code chunk number 3: HardyWeinberg.Rnw:756-757 (eval = FALSE)
###################################################
## vignette("HardyWeinberg")


###################################################
### code chunk number 4: HardyWeinberg.Rnw:771-774
###################################################
library("HardyWeinberg")
x <- c(MM = 298, MN = 489, NN = 213)
HW.test <- HWChisq(x, verbose = TRUE)


###################################################
### code chunk number 5: HardyWeinberg.Rnw:795-796
###################################################
HW.test <- HWChisq(x, cc = 0, verbose = TRUE)


###################################################
### code chunk number 6: HardyWeinberg.Rnw:804-805
###################################################
HW.lrtest <- HWLratio(x, verbose = TRUE)


###################################################
### code chunk number 7: HardyWeinberg.Rnw:813-814
###################################################
HW.exacttest <- HWExact(x, verbose = TRUE)


###################################################
### code chunk number 8: HardyWeinberg.Rnw:833-835
###################################################
set.seed(123)
HW.permutationtest <- HWPerm(x, verbose = TRUE)


###################################################
### code chunk number 9: HardyWeinberg.Rnw:849-851
###################################################
x <- c(MN = 489, NN = 213, MM = 298)
HW.test <- HWChisq(x, verbose = TRUE)


###################################################
### code chunk number 10: HardyWeinberg.Rnw:869-870
###################################################
HW.results <- HWAlltests(x, verbose = TRUE, include.permutation.test = TRUE)


###################################################
### code chunk number 11: HardyWeinberg.Rnw:881-883
###################################################
data(Markers)
Markers[1:12,]


###################################################
### code chunk number 12: HardyWeinberg.Rnw:896-900
###################################################
Xt <- table(Markers[,1])
Xv <- as.vector(Xt)
names(Xv) <- names(Xt)
HW.test <- HWChisq(Xv,cc=0,verbose=TRUE)


###################################################
### code chunk number 13: HardyWeinberg.Rnw:912-914
###################################################
set.seed(123)
Results <- HWMissing(Markers[,1], m = 50, method = "sample", verbose=TRUE)


###################################################
### code chunk number 14: HardyWeinberg.Rnw:932-934
###################################################
set.seed(123)
Results <- HWMissing(Markers[, 1:5], m = 50, verbose = TRUE)


###################################################
### code chunk number 15: HardyWeinberg.Rnw:943-945
###################################################
set.seed(123)
Results <- HWMissing(Markers[, 1:5], m = 50, statistic = "exact", verbose = TRUE)


###################################################
### code chunk number 16: HardyWeinberg.Rnw:957-959
###################################################
SNP1 <- c(A=399,B=205,AA=230,AB=314,BB=107) 
HWChisq(SNP1,cc=0,x.linked=TRUE,verbose=TRUE)


###################################################
### code chunk number 17: HardyWeinberg.Rnw:964-965
###################################################
HWChisq(SNP1[3:5],cc=0)


###################################################
### code chunk number 18: HardyWeinberg.Rnw:973-974
###################################################
HWExact(SNP1,x.linked=TRUE)


###################################################
### code chunk number 19: HardyWeinberg.Rnw:979-980
###################################################
HWExact(SNP1,x.linked=TRUE,pvaluetype="midp")


###################################################
### code chunk number 20: HardyWeinberg.Rnw:986-987
###################################################
HWExact(SNP1[3:5])


###################################################
### code chunk number 21: HardyWeinberg.Rnw:992-993
###################################################
HWPerm(SNP1,x.linked=TRUE)


###################################################
### code chunk number 22: HardyWeinberg.Rnw:998-999
###################################################
HWLratio(SNP1,x.linked=TRUE)


###################################################
### code chunk number 23: HardyWeinberg.Rnw:1004-1005
###################################################
HWAlltests(SNP1,x.linked=TRUE,include.permutation.test=TRUE)


###################################################
### code chunk number 24: HardyWeinberg.Rnw:1010-1011
###################################################
AFtest(SNP1)


###################################################
### code chunk number 25: HardyWeinberg.Rnw:1021-1022
###################################################
HWPosterior(SNP1)


###################################################
### code chunk number 26: HardyWeinberg.Rnw:1047-1056
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
### code chunk number 27: HardyWeinberg.Rnw:1129-1152
###################################################
set.seed(123)
n <- 100
m <- 100
X1 <- HWData(m, n, p = rep(0.5, m))
X2 <- HWData(m, n)
X3 <- HWData(m, n, p = rep(0.5, m), f = rep(0.5, m))
X4 <- HWData(m, n, f = rep(0.5, m))
X5 <- HWData(m, n, p = rep(c(0.2, 0.4, 0.6, 0.8), 25), pfixed = TRUE)
X6 <- HWData(m, n, exactequilibrium = TRUE)
opar <- par(mfrow = c(3, 2),mar = c(1, 0, 3, 0) + 0.1)
par(mfg = c(1, 1))
HWTernaryPlot(X1, main = "(a)", vbounds = FALSE)
par(mfg = c(1, 2))
HWTernaryPlot(X2, main = "(b)", vbounds = FALSE)
par(mfg = c(2, 1))
HWTernaryPlot(X3, main = "(c)", vbounds = FALSE)
par(mfg = c(2, 2))
HWTernaryPlot(X4, main = "(d)", vbounds = FALSE)
par(mfg = c(3, 1))
HWTernaryPlot(X5, main = "(e)", vbounds = FALSE)
par(mfg = c(3, 2))
HWTernaryPlot(X6, main = "(f)", vbounds = FALSE)
par(opar)


###################################################
### code chunk number 28: HardyWeinberg.Rnw:1177-1180 (eval = FALSE)
###################################################
## data("HapMapCHBChr1", package = "HardyWeinberg")
## HWTernaryPlot(HapMapCHBChr1, region = 1, vbounds = FALSE)
## HWTernaryPlot(HapMapCHBChr1, region = 7, vbounds = FALSE)


###################################################
### code chunk number 29: HardyWeinberg.Rnw:1213-1220 (eval = FALSE)
###################################################
## set.seed(123)
## data("HapMapCHBChr1", package = "HardyWeinberg")
## HWQqplot(HapMapCHBChr1)
## dev.off()
## set.seed(123)
## SimulatedData <- HWData(nm = 225, n = 84, p = af(HapMapCHBChr1))$Xt
## HWQqplot(SimulatedData)


