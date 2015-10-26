rStrauss.sp2 <- function(beta, gamma, R=0, p=0.5, m=100, win=sphwin()) {
stopifnot(inherits(win, "sphwin") && beta > 0 && gamma >=0 && gamma <= 1 && R >= 0 && p >=0 && p <= 1 && m > 0)
if(gamma==0) {print("Warning: Since gamma=0, simulated pattern from a Hard-Core process, not a Strauss process.")}
X <- rpoispp.sp2(beta, win=win, as.sp=FALSE)
for (i in 1:m) {
n <- nrow(X)
prop <- runif(1,0,1)
if(prop <= p) {X <- birth.sp2(X=X, beta=beta, gamma=gamma, R=R, p=p, n=n, win=win)} else {
X <- death.sp2(X=X, beta=beta, gamma=gamma, R=R, p=p, n=n, win=win)
}
}
output <- sp2(X=X, win=win)
output
}
