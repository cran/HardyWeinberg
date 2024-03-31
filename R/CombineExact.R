CombineExact <-function (Xmat, alternative="less")
{
    m <- nrow(Xmat)
    pvecless <- HWExactStats(Xmat,plinkcode=FALSE,pvaluetype="midp",alternative="less")
#    pvecgrea <- HWExactMat(Xmat,pvaluetype="midp",alternative="greater")$pvalvec
#    pvecgrea <- 1 - pvecless
    mipless <- mipvalue(pvecless)
#    mipgrea <- mipvalue(pvecgrea)
    mipgrea <- 1 - mipless
    return(list(mipless=mipless,mipgrea=mipgrea))
}
