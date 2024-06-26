HWTernaryPlot <- function(X, n=NA, addmarkers=TRUE, newframe=TRUE, hwcurve=TRUE, vbounds=FALSE, mafbounds=FALSE, mafvalue=0.05, axis=0, region=1, vertexlab=colnames(X),
                          alpha = 0.05, vertex.cex = 1, pch = 19, cc = 0.5, markercol = "black", markerbgcol= "black", cex=0.75, axislab ="", verbose=FALSE,
			  markerlab=NULL, markerpos=NULL, mcex=1, connect = FALSE, curvecols=rep("black",5), signifcolour = TRUE, patternsigsymbol = 19,
			  curtyp = "solid", ssf = "max", pvaluetype = "selome", grid = FALSE, gridlabels = TRUE, patternramp = FALSE, axisticklabels = FALSE, ...) {
#
# plot a ternary diagram that represents all rows of X as points.
# Acceptance regions for various tests for HWE can be added.
#
  if(is.vector(X)) {
    if(length(X)!=3) {
      stop("X must have three elements")
    } else {
       X <- matrix(X,ncol=3,dimnames=list(c("1"),names(X)))
    }
  }
  X <- as.matrix(X)
  nr <- nrow(X)
  nc <- ncol(X)
  if(any(X<0)) stop("X must be non-negative")
  if(nc != 3) stop("X must have three columns")
  if(is.na(n)) {
    if((sum(apply(X,1,sum))==nr)) { # data are compositions
      stop("argument n (the sample size) should be supplied")
    } else { # raw counts
      ssf <- match.fun(ssf)
      vsums <- as.matrix(apply(X,1,sum),ncol=1)
      n <- apply(vsums,2,ssf)
      Xr <- X
      if (nrow(X) == 1) { 
        Xcom <- X/sum(X)
      } else {
        Xcom <- HWClo(X)
      }
    }
  } else {
    if((sum(apply(X,1,sum))==nr)) {  # data are compositions
      Xr <- round(n*X)
      Xcom <- X
    } else { # raw counts
      Xr <- X
      if (nrow(X) == 1) {
        Xcom <- X/sum(X)
      } else {
        Xcom <- HWClo(X)
      }
    }
  }
  chiquant <- qchisq(1-alpha,1)
  r <- sqrt(chiquant/n)
  k <- cc/n
  M <- matrix(c(-1/sqrt(3),0,0,1,1/sqrt(3),0),ncol=2,byrow=T)
  #
  # coordinates for grid
  #
  gridstyle <- "dashed"
  B1 <- rbind(c(0.8,0.2,0.0),
              c(0.6,0.4,0.0),
              c(0.4,0.6,0.0),
              c(0.2,0.8,0.0))
  E1 <- rbind(c(0.0,0.2,0.8),
              c(0.0,0.4,0.6),
              c(0.0,0.6,0.4),
              c(0.0,0.8,0.2))
  B2 <- B1[,c(3,1,2)]
  E2 <- E1[,c(3,1,2)]
  B3 <- B1[,c(2,3,1)]
  E3 <- E1[,c(2,3,1)]
  B <- rbind(B1,B2,B3)
  E <- rbind(E1,E2,E3)
  Bc <- B%*%M
  Ec <- E%*%M
  nsignif <- NA
  markerq <- (Xcom[,2]+2*Xcom[,3])/2
  if(newframe) {
    opar <- par(pty="m",xpd=TRUE)
    on.exit(par(opar))
    plot(M[,1],M[,2], type="n", axes=FALSE, xlab="", ylab="", pch=19, asp=1, cex.main=2, ... )
    polygon(M)
    eps <- 0.04 * vertex.cex
    Mlab <- M + matrix(c(-eps,0,0,eps,eps,0),ncol=2,byrow=T)
    text(Mlab[,1],Mlab[,2], vertexlab, cex=vertex.cex)
    text(0,-0.1,axislab,cex=vertex.cex)
  }
  if (axis==0) {
    ;
  } else if (axis==1) {
    AXA <- rbind(c(0,0.5,0.5),c(1,0,0))
    AXA <- AXA%*%M
    lines(AXA[,1],AXA[,2],...)      
  } else if (axis==2) {
    AXAB <- rbind(c(0.5,0,0.5),c(0,1,0))
    AXAB <- AXAB%*%M
    lines(AXAB[,1],AXAB[,2],...)
  } else if (axis==3) {
    AXB <- rbind(c(0.5,0.5,0),c(0,0,1))
    AXB <- AXB%*%M
    lines(AXB[,1],AXB[,2],...)
  } else if (axis==4) {
    l <- 2/sqrt(3) # total length of base
    steps <- l/10
    ss <- seq(-1/sqrt(3),1/sqrt(3),by=steps)
    fr <- cbind(ss,rep(0,11))
    to <- cbind(ss,rep(-0.05,11))
    segments(fr[,1],fr[,2],to[,1],to[,2])
    if(axisticklabels) {
      for(i in 0:10) {
        text(to[i+1,1],to[i+1,2],toString(round(i/10,digits=1)),pos=1,cex=0.75)
      }
    }
  } else {
    cat("Unknown value for parameter axis\n")
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
    
  D <- 0.5*(Xcom[,2] - 2*(1-markerq)*markerq)
  Dpos <- sum(Xcom[,2] > 2*(1-markerq)*markerq)
  Dneg <- sum(Xcom[,2] < 2*(1-markerq)*markerq)
  Dzer <- sum(Xcom[,2] == 2*(1-markerq)*markerq)
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
  if (grid) {
    for(i in 1:nrow(Bc)) {
      segments(Bc[i,1],Bc[i,2],Ec[i,1],Ec[i,2],lty=gridstyle)
    }

    if(gridlabels) {
    
      text(Ec[1,1],Ec[1,2],"0.2",pos=4,cex=0.75)
      text(Ec[2,1],Ec[2,2],"0.4",pos=4,cex=0.75)
      text(Ec[3,1],Ec[3,2],"0.6",pos=4,cex=0.75)
      text(Ec[4,1],Ec[4,2],"0.8",pos=4,cex=0.75)
    
      text(Ec[5,1],Ec[5,2],"0.2",pos=1,cex=0.75,srt=-120,offset=1)
      text(Ec[6,1],Ec[6,2],"0.4",pos=1,cex=0.75,srt=-120,offset=1)
      text(Ec[7,1],Ec[7,2],"0.6",pos=1,cex=0.75,srt=-120,offset=1)
      text(Ec[8,1],Ec[8,2],"0.8",pos=1,cex=0.75,srt=-120,offset=1)
    
      text(Ec[9,1],Ec[9,2],"0.2",pos=3,cex=0.75,srt=+120,offset=1)
      text(Ec[10,1],Ec[10,2],"0.4",pos=3,cex=0.75,srt=+120,offset=1)
      text(Ec[11,1],Ec[11,2],"0.6",pos=3,cex=0.75,srt=+120,offset=1)
      text(Ec[12,1],Ec[12,2],"0.8",pos=3,cex=0.75,srt=+120,offset=1)
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

  if(region==0) {
    ;
  } else if(region==1) { # simple hw ci
    HWChisqUpperl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2],curtyp=curtyp)
    HWChisqLowerl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2],curtyp=curtyp)
  } else if(region==2) { # all curves for hw with cc
    # lower curve for D<0

    DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)

    # upper curve for D<0

    DnegUL <- DnegUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)

    # lower curve for D>0

    DposLL <- DposLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)
      
    # upper curve for D>0

    DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)

    if(verbose) {
      cat("D<0 LL",round(DnegLL,digits=2),"\n")
      cat("D<0 UL",round(DnegUL,digits=2),"\n")
      cat("D>0 LL",round(DposLL,digits=2),"\n")
      cat("D>0 UL",round(DposUL,digits=2),"\n")
    }
  } else if(region==3) { # only limits for D>0, chisq with cc
    # lower curve for D>0

    DposLL <- DposLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)
      
    # upper curve for D>0

    DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)

  } else if(region==4) { # only limits for D<0, chisq with cc
    # lower curve for D<0

    DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)

    # upper curve for D<0

    DnegUL <- DnegUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)
  } else if(region==5) { # all limits
    HWChisqUpperl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])
    HWChisqLowerl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])

    # lower curve for D<0

    DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)

    # upper curve for D<0

    DnegUL <- DnegUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)

    # lower curve for D>0

    DposLL <- DposLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)
      
    # upper curve for D>0

    DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)
         
  } else if(region==6) { # lower for D<0, upper for D>0

    #         HWChisqUpperl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])
    #         HWChisqLowerl(r,verbose=FALSE,cex=cex,curvecol=curvecols[2])

    # lower curve for D<0

    DnegLL <- DnegLow(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[4],curtyp=curtyp)

    # upper curve for D>0

    DposUL <- DposUp(r,k,n,cc,chiquant,verbose=verbose,cex=cex,curcol=curvecols[3],curtyp=curtyp)

  } else if(region==7) { # For Haldane's Exact test
    Crit <- CritSam(n,alphalimit=alpha,pvaluetype=pvaluetype)$Xn
    Critcar <- Crit%*%M        # cartesian coordinates
    points(Critcar[,1],Critcar[,2],pch=pch,col=curvecols[5],cex=cex,type="l",lty=curtyp)

    Crit <- CritSam(n,Dpos=FALSE,alphalimit=alpha,pvaluetype=pvaluetype)$Xn
    Critcar <- Crit%*%M        # cartesian coordinates
    points(Critcar[,1],Critcar[,2],pch=pch,col=curvecols[5],cex=cex,type="l",lty=curtyp)
  } else {
    cat("Unknown value for parameter region\n")
  }
  if(addmarkers) {
    Xc <- Xcom%*%M # cartesian coordinates
    if (patternramp) {
      Xu <- UniqueGenotypeCounts(X,verbose=FALSE)
      w <- Xu[,4]/sum(Xu[,4])
      wr <- w/max(w)
      colvec <- MakeColourVector(wr)
      ii <- order(Xu[,4],decreasing=TRUE)
      Xu <- Xu[ii,]
      colvec <- colvec[ii]
      Xcom <- HWClo(as.matrix(Xu[,1:3]))
      Xc <- Xcom%*%M # cartesian coordinates
      # use a different plotting symbol for significant markers if required
      if(patternsigsymbol != pch) {
        npatterns   <- nrow(Xc)
        markersym   <- rep(pch,npatterns)   
        if (region == 1) {
          chi.stats   <- numeric(npatterns)           
          chisq.crit  <- qchisq(1-alpha,1)
          chisq.stats <- HW.chi.mat(Xu) 
          markersym[chisq.stats > chisq.crit] <- patternsigsymbol
        } else if (region == 2) {
	  options(warn=-1)
          pvals <- numeric(npatterns)
          for (i in 1:npatterns) {
            pvals[i] <- HWChisq(Xu[i,],cc=0.5,verbose=FALSE)$pval
          }
          markersym[pvals<alpha] <- patternsigsymbol
	  options(warn=0)
        } else if (region == 7) {
          pvals <- numeric(npatterns)
          for (i in 1:npatterns) {
            pvals[i] <- HWExact(Xu[i,],alternative="two.sided",verbose=FALSE,pvaluetype)$pval
          }
          markersym[pvals<alpha] <- patternsigsymbol
        } else {
          ;
        }
        points(Xc[,1],Xc[,2],pch=markersym,bg=markerbgcol,col=colvec,cex=cex)  
      } else {
        points(Xc[,1],Xc[,2],pch=pch,bg=markerbgcol,col=colvec,cex=cex)  
      }
      # add the colour ramp
      colpal<- colorRampPalette(c("green","yellow", "red"))(n = 1000)
      colorlegend(colpal,zlim=c(0,max(w)),posx=c(0.85,0.88),
            posy=c(0.55,0.95),zval=c(0.00,max(w)),
            digit=3,cex=0.50)
    } else {
      if (signifcolour==TRUE) { 
        if (region == 1) {
          chi.stats <- numeric(nr)           
          chisq.crit <- qchisq(1-alpha,1)
          chisq.stats <- HW.chi.mat(Xr) 
          markerbgcol <- rep("green",nr)             
          markerbgcol[chisq.stats > chisq.crit] <- "red" # look if chisquare is too large.
          markercol <- rep("green",nr)             
          markercol[chisq.stats > chisq.crit] <- "red"
          nsignif <- sum(chisq.stats > chisq.crit)
        } else if (region == 2) {
	  options(warn=-1)
          pvals <- numeric(nr)
          for (i in 1:nr)
            pvals[i] <- HWChisq(Xr[i,],cc=0.5,verbose=FALSE)$pval
          markerbgcol <- rep("green",nr)             
          markerbgcol[pvals<alpha] <- "red"
          markercol <- rep("green",nr)             
          markercol[pvals<alpha] <- "red"
          nsignif <- sum(pvals<alpha)
	  options(warn=0)
        } else if (region == 7) {
          pvals <- numeric(nr)
          for (i in 1:nr) {
            x <- Xr[i,]
            pvals[i] <- HWExact(Xr[i,],alternative="two.sided",verbose=FALSE,pvaluetype)$pval
            markerbgcol <- rep("green",nr)             
            markerbgcol[pvals<alpha] <- "red"
            markercol <- rep("green",nr)             
            markercol[pvals<alpha] <- "red"
            nsignif <- sum(pvals<alpha)
          }
        } else {
          ;
        }
      } # end if signifcolour
      if (connect) {
        points(Xc[,1],Xc[,2],pch=pch,col=curvecols[1],cex=cex,type="l",lty=curtyp)
        points(Xc[,1],Xc[,2],pch=pch,bg=markerbgcol,col=markercol,cex=cex)
      } else {
        points(Xc[,1],Xc[,2],pch=pch,bg=markerbgcol,col=markercol,cex=cex)
	text(Xc[,1],Xc[,2],markerlab,cex=mcex,pos=markerpos)
      }
    } # else patternramp
  } # if addmarkers
  results <- list(minp=minp,maxp=maxp,inrange=inrange,percinrange=percinrange,nsignif=nsignif)
}

