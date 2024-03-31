MakeCounts <- function (X, alleles=NULL, pos1 = 1, pos2 = 3, coding = c(AA = 0, AB = 1, BB = 2), sep = "") 
{
    if (is.vector(X)) {
        X <- as.matrix(X, ncol = 1)
    }
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    Y <- matrix(NA, nrow = p, ncol = 4)
    if (is.null(colnames(X))) {
        rownames(Y) <- paste("M", 1:p, sep = "")
    }
    else {
        rownames(Y) <- colnames(X)
    }
    colnames(Y) <- c("AA", "AB", "BB", "NA")
    if (!is.numeric(X)) {
      if(is.null(alleles)) { # extract alleles from the data
        als <- character(p)
        for(j in 1:p) {
          thisone <- alleles(X[,j],fromlabels = FALSE)
          if(length(thisone) == 2) { # biallelic
            als[j] <- paste(thisone[1],thisone[2],sep="/")
          } else if(length(thisone) == 1) { # monomorphic
            als[j] <- paste(thisone,"$",sep="/")
          } else {
            stop(paste("Column",toString(j),
                       "is a multi-allelic marker, recode as bi-allelic.\n"))
          }
        }
        alleles <- als
      }
      for (j in 1:p) {
          snp <- X[, j]
          al1 <- substr(alleles[j], pos1, pos1)
          al2 <- substr(alleles[j], pos2, pos2)
          homAA <- paste(al1, al1, sep = sep)
          homBB <- paste(al2, al2, sep = sep)
          hetAB <- paste(al1, al2, sep = sep)
          hetBA <- paste(al2, al1, sep = sep)
          nAA <- sum(snp == homAA, na.rm = TRUE)
          nBB <- sum(snp == homBB, na.rm = TRUE)
          nAB <- sum(snp == hetAB, na.rm = TRUE) + sum(snp == 
                hetBA, na.rm = TRUE)
          nNA <- sum(is.na(snp))
          tot <- nAA + nAB + nBB + nNA
          if (tot != n) {
            cat(j, "\n")
            cat("nAA =", nAA, "nAB =", nAB, "nBB =", nBB, 
                  "nNA =", nNA, "n = ", n, "\n")
            stop("genotypes and missings do not sum n")
          }
          Y[j, ] <- c(nAA, nAB, nBB, nNA)
        } # end for
    }
    else {
        for (j in 1:p) {
            snp <- X[, j]
            homAA <- coding[1]
            homBB <- coding[3]
            hetAB <- coding[2]
            nAA <- sum(snp == homAA, na.rm = TRUE)
            nBB <- sum(snp == homBB, na.rm = TRUE)
            nAB <- sum(snp == hetAB, na.rm = TRUE)
            nNA <- sum(is.na(snp))
            tot <- nAA + nAB + nBB + nNA
            if (tot != n) {
                cat(j, "\n")
                cat("nAA =", nAA, "nAB =", nAB, "nBB =", nBB, 
                  "nNA =", nNA, "n = ", n, "\n")
                stop("genotypes and missings do not sum n")
            }
            Y[j, ] <- c(nAA, nAB, nBB, nNA)
        }
    }
    return(Y)
}
