rHardcore.sp2 <- function(beta, R=0, p=0.5, m=100, win=sphwin(), proper=TRUE) {
  stopifnot(inherits(win, "sphwin"))
  X <- rpoispp.sp2(beta, as.sp=FALSE)
  n <- nrow(X)
  if(proper) {
    while(sort(gcdist(X,X))[n+1] <= R) {
      prop <- runif(1,0,1)
      if(prop <= p) {
        X1 <- birth.sp2(X=X, beta=beta, gamma=0, R=R, p=p, n=n, win=win)
      } else {
        X1 <- death.sp2(X=X, beta=beta, gamma=0, R=R, p=p, n=n, win=win)
      }
      X <- X1
      n <- nrow(X)
    }
  } else {
    for(i in 1:m) {
      prop <- runif(1,0,1)
      if(prop <= p) {
        X1 <- birth.sp2(X=X, beta=beta, gamma=0, R=R, p=p, n=n, win=win)
      } else {
        X1 <- death.sp2(X=X, beta=beta, gamma=0, R=R, p=p, n=n, win=win)
      }
      X <- X1
      n <- nrow(X)
    }
  }
  output <- sp2(X=X, win=win)
  output
}

