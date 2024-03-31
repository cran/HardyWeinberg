toTriangular <- function(x) 
{
    x   <- unlist(x)
    k   <- n.alleles(x)
    als <- alleles(x)
    als.num <- 1:k
    names(als.num) <- als
    A1 <- substr(names(x),1,1)
    A2 <- substr(names(x),2,2)
    i1 <- rep(A1,x)
    i2 <- rep(A2,x)
    GTc <- AllelesToTriangular(i1,i2)
    nallelesfound <- nrow(GTc)
    if(k!=nallelesfound) {
      stop("Specified alleles are not all found.\n")  
    } else {
      colnames(GTc) <- als
      rownames(GTc) <- als   
    }
    return(GTc)
}
