HWAIC <- function(x, y, tracing = 0) {
  z <- x+y
  pam <- af(x)
  paf <- af(y)
  pa  <- af(z)
  fm <- HWf(x)
  ff <- HWf(y)
  f  <- HWf(z)
  l1 <- loglik.1(pa,z) # EoAF and HWP (f's zero) (S1)
  l2 <- loglik.2(pa,z,f) # EoAF and EoIC (S2)
  l3 <- loglik.3(pa,f,x,y,tracing=tracing) # EAF only (S3)
  l4 <- loglik.4(pa,x,y,pam,paf) # HWP for both sexes (S4)
  l5 <- loglik.5(x,y,tracing=tracing) #no EAF but with EIC (S5)
  l6 <- loglik.6(x,y,pam,paf,fm,ff)
  ll <- c(l1[1],l2[1],l3[1],l4[1],l5[1],l6[1])
  np <- c(1,2,3,2,3,4)
  aic <- 2*np-2*ll
  names(aic) <- c("EAF & HWP","EAF & EIC","EAF","HWP","EIC","NONE")
  return(aic)
}
