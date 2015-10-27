birth.sp2 <- function(X, beta, gamma, R, p, n, rad=1, win) {
rad <- win$rad
win1 <- win$win
w <- runif.sphwin(n=1,win=win, as.sp=FALSE)
x <- convert3(w)
y <- convert3(X)
xs <- ifelse(gcdist(x,y) <= R, 1, 0)
birthcheck <- beta*(gamma^(sum(xs)))*((1-p)/(n+1))*(4*pi*rad^2)/p
probbirth <- min(1, birthcheck)
y <- runif(1,0,1)
if(y >= probbirth) {Xfinal <- X} else {Xfinal <- rbind(X, w)}
Xfinal
}
