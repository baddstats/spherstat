rMaternI.sp2 <- function(kappa, r, win=sphwin(type="sphere")) {
X1 <- rpoispp.sp2(lambda=kappa, win=win)
nrX <- nrow(X1$X)
X2 <- runif.sphwin(n=1, win=win)
for(i in 1:(nrX-1)) {
Xtest <- runif.sphwin(n=1, win=win)
if(min(nncrosssph(Xtest, X2)) > r) {X2$X <- rbind(X2$X, Xtest$X)}
}
X2
}