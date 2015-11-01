birth.sp2 <- function(X, beta, gamma, R, p, n, rad=1, win) {
  rad <- win$rad
  w <- runif.sphwin(n=1,win=win, as.sp=FALSE)
  x <- convert3(w)
  y <- convert3(X)
  xs <- (gcdist(x,y) <= R)
  birthcheck <- beta*(gamma^(sum(xs)))*((1-p)/(n+1))*(4*pi*rad^2)/p
  probbirth <- min(1, birthcheck)
  y <- runif(1,0,1)
  if(y >= probbirth) {
    ## proposal rejected
    return(X)
  }
  ## proposal accepted
  return(rbind(X, w))
}
