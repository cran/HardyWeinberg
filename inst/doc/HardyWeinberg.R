### R code from vignette source 'HardyWeinberg.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: HardyWeinberg.Rnw:743-744
###################################################
options(prompt = "R> ", continue = "+ ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: HardyWeinberg.Rnw:747-749 (eval = FALSE)
###################################################
## install.packages("HardyWeinberg")
## library("HardyWeinberg")


###################################################
### code chunk number 3: HardyWeinberg.Rnw:760-761 (eval = FALSE)
###################################################
## vignette("HardyWeinberg")


###################################################
### code chunk number 4: HardyWeinberg.Rnw:775-778
###################################################
library("HardyWeinberg")
x <- c(MM = 298, MN = 489, NN = 213)
HW.test <- HWChisq(x, verbose = TRUE)


###################################################
### code chunk number 5: HardyWeinberg.Rnw:799-800
###################################################
HW.test <- HWChisq(x, cc = 0, verbose = TRUE)


###################################################
### code chunk number 6: HardyWeinberg.Rnw:808-809
###################################################
HW.lrtest <- HWLratio(x, verbose = TRUE)


###################################################
### code chunk number 7: HardyWeinberg.Rnw:817-818
###################################################
HW.exacttest <- HWExact(x, verbose = TRUE)


###################################################
### code chunk number 8: HardyWeinberg.Rnw:837-839
###################################################
set.seed(123)
HW.permutationtest <- HWPerm(x, verbose = TRUE)


###################################################
### code chunk number 9: HardyWeinberg.Rnw:853-855
###################################################
x <- c(MN = 489, NN = 213, MM = 298)
HW.test <- HWChisq(x, verbose = TRUE)


###################################################
### code chunk number 10: HardyWeinberg.Rnw:873-874
###################################################
HW.results <- HWAlltests(x, verbose = TRUE, include.permutation.test = TRUE)


###################################################
### code chunk number 11: HardyWeinberg.Rnw:885-887
###################################################
data(Markers)
Markers[1:12,]


###################################################
### code chunk number 12: HardyWeinberg.Rnw:900-904
###################################################
Xt <- table(Markers[,1])
Xv <- as.vector(Xt)
names(Xv) <- names(Xt)
HW.test <- HWChisq(Xv,cc=0,verbose=TRUE)


###################################################
### code chunk number 13: HardyWeinberg.Rnw:916-918
###################################################
set.seed(123)
Results <- HWMissing(Markers[,1], m = 50, method = "sample", verbose=TRUE)


###################################################
### code chunk number 14: HardyWeinberg.Rnw:936-938
###################################################
set.seed(123)
Results <- HWMissing(Markers[, 1:5], m = 50, verbose = TRUE)


###################################################
### code chunk number 15: HardyWeinberg.Rnw:947-949
###################################################
set.seed(123)
Results <- HWMissing(Markers[, 1:5], m = 50, statistic = "exact", verbose = TRUE)


###################################################
### code chunk number 16: HardyWeinberg.Rnw:957-959
###################################################
data(JPTsnps)
Results <- HWPosterior(JPTsnps[1,],x.linked=FALSE,precision=0.05)


###################################################
### code chunk number 17: HardyWeinberg.Rnw:966-969
###################################################
data(JPTsnps)
AICs <- HWAIC(JPTsnps[1,1:3],JPTsnps[1,4:6])
AICs


###################################################
### code chunk number 18: HardyWeinberg.Rnw:981-983
###################################################
SNP1 <- c(A=399,B=205,AA=230,AB=314,BB=107) 
HWChisq(SNP1,cc=0,x.linked=TRUE,verbose=TRUE)


###################################################
### code chunk number 19: HardyWeinberg.Rnw:988-989
###################################################
HWChisq(SNP1[3:5],cc=0)


###################################################
### code chunk number 20: HardyWeinberg.Rnw:997-998
###################################################
HWExact(SNP1,x.linked=TRUE)


###################################################
### code chunk number 21: HardyWeinberg.Rnw:1003-1004
###################################################
HWExact(SNP1,x.linked=TRUE,pvaluetype="midp")


###################################################
### code chunk number 22: HardyWeinberg.Rnw:1010-1011
###################################################
HWExact(SNP1[3:5])


###################################################
### code chunk number 23: HardyWeinberg.Rnw:1016-1017
###################################################
HWPerm(SNP1,x.linked=TRUE)


###################################################
### code chunk number 24: HardyWeinberg.Rnw:1022-1023
###################################################
HWLratio(SNP1,x.linked=TRUE)


###################################################
### code chunk number 25: HardyWeinberg.Rnw:1028-1029
###################################################
HWAlltests(SNP1,x.linked=TRUE,include.permutation.test=TRUE)


###################################################
### code chunk number 26: HardyWeinberg.Rnw:1034-1035
###################################################
AFtest(SNP1)


###################################################
### code chunk number 27: HardyWeinberg.Rnw:1045-1046
###################################################
HWPosterior(SNP1,x.linked=TRUE)


###################################################
### code chunk number 28: HardyWeinberg.Rnw:1071-1080
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
### code chunk number 29: HardyWeinberg.Rnw:1153-1176
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
### code chunk number 30: HardyWeinberg.Rnw:1201-1204 (eval = FALSE)
###################################################
## data("HapMapCHBChr1", package = "HardyWeinberg")
## HWTernaryPlot(HapMapCHBChr1, region = 1, vbounds = FALSE)
## HWTernaryPlot(HapMapCHBChr1, region = 7, vbounds = FALSE)


###################################################
### code chunk number 31: HardyWeinberg.Rnw:1237-1244 (eval = FALSE)
###################################################
## set.seed(123)
## data("HapMapCHBChr1", package = "HardyWeinberg")
## HWQqplot(HapMapCHBChr1)
## dev.off()
## set.seed(123)
## SimulatedData <- HWData(nm = 225, n = 84, p = af(HapMapCHBChr1))$Xt
## HWQqplot(SimulatedData)


