in.W.poly <-
function(points, win=sphwin(type="sphere")) {
stopifnot(inherits(points, c("sp2", "sp3", "matrix")))
if(inherits(points, "matrix")) {stopifnot(inherits(win, "sphwin"))} else {
win <- points$win
points <- points$X
stopifnot(inherits(win, "sphwin"))
}
rad <- win$rad
n <- nrow(points)
lp <- nrow(win$param)
p3 <- convert3(win$param)
if(ncol(points) !=3) {points <- convert3(points)}
polyangs1 <- p3[1:(lp-1),]
polyangs2 <- p3[2:lp,]
polyangs3 <- rbind(p3[3:lp, ], p3[2,])
gc <- diag(gcdist(x=polyangs1, y=polyangs2, rad=rad))
gc1 <- gc[1:(lp-1)]
gc2 <- c(gc[2:(lp-1)], gc[1])
sphang <- sph.angles(win)
gc.pv1 <- gcdist(x=points, y=polyangs1)
gc.pv2 <- gcdist(x=points, y=polyangs2)
gc.pv3 <- gcdist(x=points, y=polyangs3)
			sc1 <- sphcos(d1=gc1, d2=gc.pv2, d3=gc.pv1, theta=NULL, rad=rad)
			sc2 <- sphcos(d1=gc2, d2=gc.pv2, d3=gc.pv3, theta=NULL, rad=rad)
test <- abs(sc1+sc2-matrix(rep(sphang, n), nrow=n, ncol=lp-1, byrow=TRUE))
test <- ifelse(is.na(test), 0, test)
test1 <- rowSums(sround(test)<10^-7)
isinW <- test1==lp-1
isinW
}