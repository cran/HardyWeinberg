shannon <- function(x) {
  S <- length(x)
  N <- sum(x)
  x <- x/N
  x <- x[x!=0]
  Hp <- -sum(x*log(x))
  VHp <- (sum(x*log(x)^2) - Hp^2)/N + (S-1)/(2*N^2)
  return(list(Hp=Hp,VHp=VHp))
}
