`HWTernaryPlot` <-
function(X, n, addmarkers=TRUE, newframe=TRUE, hwcurve=TRUE, vbounds=TRUE, mafbounds=FALSE, mafvalue=0.05, axis=0, region=1, vertexlab=colnames(X), alpha = 0.05, vertex.cex = 1, pch = 19, cc = 0.5, markercol = "black", cex=0.75, axislab ="", verbose=FALSE, markerlab=NULL, mcex=1, connect = FALSE, curvecols=rep("black",5) , signifcolour=FALSE, ...)
  {
    X <- as.matrix(X)
    nr <- nrow(X)
    nc <- ncol(X)
    chiquant <- qchisq(1-alpha,1)
    r <- sqrt(chiquant/n)
    k <- cc/n
    M <- matrix(c(-1/sqrt(3),0,0,1,1/sqrt(3),0),ncol=2,byrow=T)
    nsignif <- NA
    
    if(any(X) < 0) stop("X must be non-negative")
    if(nc != 3) stop("X must have three columns")
    if(!(sum(apply(X,1,sum))==nr)) stop("rows of X must sum to 1")
    markerq <- (X[,2]+2*X[,3])/2

      if(newframe) {
       oldpty <- par("pty")
       on.exit(par(pty=oldpty))
       par(pty="m")
       par(xpd=TRUE) # allow for labels...

       plot(M[,1],M[,2], type="n", axes=FALSE, xlab="", ylab="", pch=19, asp=1, cex.main=2, ... )
       polygon(M)

       if (axis==1) {
          AXA <- rbind(c(0,0.5,0.5),c(1,0,0))
          AXA <- AXA%*%M
          lines(AXA[,1],AXA[,2])      
       }

       if (axis==2) {
          AXAB <- rbind(c(0.5,0,0.5),c(0,1,0))
          AXAB <- AXAB%*%M
          lines(AXAB[,1],AXAB[,2])
       }

       if (axis==3) {
          AXB <- rbind(c(0.5,0.5,0),c(0,0,1))
          AXB <- AXB%*%M
          lines(AXB[,1],AXB[,2])
       }
       
       eps <- 0.04 * vertex.cex
       Mlab <- M + matrix(c(-eps,0,0,eps,eps,0),ncol=2,byrow=T)
       text(Mlab[,1],Mlab[,2], vertexlab, cex=vertex.cex)

       text(0,-0.1,axislab,cex=vertex.cex)
     }
       if(hwcurve) {
         p <- seq(0,1,by=0.005)
         HW <- cbind(p^2,2*p*(1-p),(1-p)^2)
         HWc <- HW%*%M
         points(HWc[,1],HWc[,2],type="l",col=curvecols[1])
       }

       minp  <- sqrt(5/n)
       minpt <- 2*(minp-0.5)/sqrt(3)
       maxp  <- 1-sqrt(5/n)
       maxpt <- 2*(maxp-0.5)/sqrt(3)

       inrange <- sum((markerq >= minp & markerq <= maxp))
       percinrange <- round(100*inrange/nr,digits=2)

       ind1 <- markerq==1
       ind0 <- markerq==0
       nfixed <- sum(ind1) + sum(ind0)
    
       D <- 0.5*(X[,2] - 2*(1-markerq)*markerq)
       Dpos <- sum(X[,2] > 2*(1-markerq)*markerq)
       Dneg <- sum(X[,2] < 2*(1-markerq)*markerq)
       Dzer <- sum(X[,2] == 2*(1-markerq)*markerq)
       Dtot <- Dpos+Dneg


#       cat("D>0:",Dpos,"D<0:",Dneg,"D=0:",Dzer,"nfix:",nfixed,"Tot:",Dtot,"mD:",mean(D),"medD:",median(D),"\n")
#       cat("D>0:",round(100*Dpos/Dtot,digits=5),"D<0:",round(100*Dneg/Dtot,digits=5),"D=0:",
#           round(100*Dzer/Dtot,digits=5),"\n")
       
       if(vbounds) {
         if (n >= 20) {
            lines(c(minpt,minpt),c(0,2*minp),lty="dashed")
            lines(c(maxpt,maxpt),c(0,2-2*maxp),lty="dashed")
          }
       }

       if(mafbounds) {
          minaf <- mafvalue
          minaft <- 2*(minaf-0.5)/sqrt(3)
          maxaf <- 0.95
          maxaft <- 2*(maxaf-0.5)/sqrt(3)
          lines(c(minaft,minaft),c(0,2*minaf),lty="dashed")
          lines(c(maxaft,maxaft),c(0,2-2*maxaf),lty="dashed")
       }    

       if(region==1) { # simple hw ci
          HWChisqUpperl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])
          HWChisqLowerl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])
       }

       if(region==2) { # all curves for hw with cc

         # lower curve for D<0

         DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp="solid")

         # upper curve for D<0

         DnegUL <- DnegUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp="solid")

         # lower curve for D>0

         DposLL <- DposLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp="solid")
      
         # upper curve for D>0

         DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp="solid")

         if(verbose) {
             cat("D<0 LL",round(DnegLL,digits=2),"\n")
             cat("D<0 UL",round(DnegUL,digits=2),"\n")
             cat("D>0 LL",round(DposLL,digits=2),"\n")
             cat("D>0 UL",round(DposUL,digits=2),"\n")
         }
       }
       if(region==3) { # only limits for D>0, chisq with cc
         # lower curve for D>0

         DposLL <- DposLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp="dotted")
      
         # upper curve for D>0

         DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp="dotted")

       }
       if(region==4) { # only limits for D<0, chisq with cc
         # lower curve for D<0

         DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp="dotted")

         # upper curve for D<0

         DnegUL <- DnegUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp="dotted")
       }

       if(region==5) { # all limits

         HWChisqUpperl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])
         HWChisqLowerl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])

         # lower curve for D<0

         DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp="dotted")

         # upper curve for D<0

         DnegUL <- DnegUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp="dotted")

         # lower curve for D>0

         DposLL <- DposLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp="dotted")
      
         # upper curve for D>0

         DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp="dotted")
         
       }

       if(region==6) { # lower for D<0, upper for D>0

         HWChisqUpperl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])
         HWChisqLowerl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])

         # lower curve for D<0

         DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp="dashed")

         # upper curve for D>0

         DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp="dashed")
         
         
       }

       if(region==7) { # For Fisher-Exact test

          Crit <- CritSam(n)$Xn
          Critcar <- Crit%*%M        # cartesian coordinates
          points(Critcar[,1],Critcar[,2],pch=pch,col=curvecols[5],cex=cex,type="l")

          Crit <- CritSam(n,Dpos=FALSE)$Xn
          Critcar <- Crit%*%M        # cartesian coordinates
          points(Critcar[,1],Critcar[,2],pch=pch,col=curvecols[5],cex=cex,type="l")

       }

       if(addmarkers) {
          Xc <- X%*%M # cartesian coordinates

          if (signifcolour==TRUE) { 

             Xa <- X*n # observed absolute frequencies
             pvals <- NULL
  
             if (region == 1) {
             
                for (i in 1:nrow(Xa))
                   pvals <- c(pvals,HWChisq(Xa[i,])$pval)
                markercol <- rep("green",nr)             
                markercol[pvals<0.05] <- "red"
                nsignif <- sum(pvals<0.05)

             }

             if (region == 2) {

                for (i in 1:nrow(Xa))
                   pvals <- c(pvals,HWChisq(Xa[i,],cc=0.5)$pval)
                markercol <- rep("green",nr)             
                markercol[pvals<0.05] <- "red"
                nsignif <- sum(pvals<0.05)

             }

             if (region == 7) {

                for (i in 1:nrow(Xa)) {
                    
                   x <- Xa[i,]
                   m <- matrix(c(x[1],x[2]/2,x[2]/2,x[3]),ncol=2)
                   out <- fisher.test(m,alternative="two.sided")
                   pvals <- c(pvals,out$p.value)
                   markercol <- rep("green",nr)             
                   markercol[pvals<0.05] <- "red"
                   nsignif <- sum(pvals<0.05)
             }

             }


          }

          if (connect)
              points(Xc[,1],Xc[,2],pch=pch,bg=markercol,col=markercol,cex=cex,type="l")
          else
              points(Xc[,1],Xc[,2],pch=pch,bg=markercol,col=markercol,cex=cex)
          text(Xc[,1],Xc[,2],markerlab,cex=mcex)
       }

       return(list(minp=minp,maxp=maxp,inrange=inrange,percinrange=percinrange,nsignif=nsignif))
  }

