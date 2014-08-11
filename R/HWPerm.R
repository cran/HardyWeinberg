HWPerm <- function (x, nperm = 17000, verbose = TRUE, FUN = Chisquare, ...) 
{
    n <- sum(x)
    nA <- 2 * x[1] + x[2]
    nB <- 2 * n - nA
    stat.obs <- FUN(x)
    pseudodist <- numeric(nperm)
    i1 <- seq(1, 2 * n, 2)
    i2 <- seq(2, 2 * n, 2)
    for (i in 1:nperm) {
        xx <- sample(c(rep("A", nA), rep("B", nB)))
        A1 <- xx[i1]
        A2 <- xx[i2]
        Geno <- paste(A1, A2, sep = "")
        Geno[Geno == "BA"] <- "AB"
        nAA <- sum(Geno == "AA")
        nAB <- sum(Geno == "AB")
        nBB <- sum(Geno == "BB")
        y <- c(AA=nAA, AB=nAB, BB=nBB)
        stat.pseudo <- FUN(y)
        pseudodist[i] <- stat.pseudo
    }
    nlarger <- sum(pseudodist >= stat.obs)
    pval <- nlarger/nperm
    if (verbose) {
        cat("Permutation test for Hardy-Weinberg equilibrium\n")
        cat(nperm, "permutations. p-value:", pval, "\n")
    }
    return(pval)
}

