rMaternII.sp2 <- function(kappa, r, win=sphwin(type="sphere")) {
X1 <- rpoispp.sp2(lambda=kappa, win=win)
nrX <- nrow(X1$X)
X2 <- runif.sphwin(n=1, win=win)
for(i in 1:(nrX-1)) {
Xtest <- runif.sphwin(n=1, win=win)
to.keep <- c()
nrX2 <- nrow(X2$X)
if(length(X2$X)==0) {cat("0 ")
X2$X <- rbind(X2$X, Xtest)} else {
for(j in 1:nrX2) {
if(gcdist(Xtest$X, matrix(X2$X[j,], nrow=1, ncol=2, byrow=TRUE)) > r) {
to.keep <- c(to.keep, j)}
}
if(length(to.keep) == nrX2) {
X2$X <- rbind(X2$X, Xtest$X)} else {
X2$X <- X2$X[to.keep,]}
}
}
X2
}

