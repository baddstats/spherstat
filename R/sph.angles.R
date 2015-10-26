sph.angles <-
function (win) {
stopifnot(inherits(win, c("sp2", "sp3", "sphwin")))
if(inherits(win, "sphwin")) {win <- win} else {win <- win$win}
rad <- win$rad
wp3 <- convert3(win$param)
lp <- nrow(wp3)
vert1 <- wp3[1:(lp-1),]
vert2 <- wp3[2:lp,]
vert3 <- rbind(wp3[3:lp,], wp3[2,])
dist1 <- diag(gcdist(x=vert1, y=vert2, rad=rad))
dist2 <- diag(gcdist(x=vert2, y=vert3, rad=rad))
dist3 <- diag(gcdist(x=vert1, y=vert3, rad=rad))
angles <- mapply(sphcos, d1=dist1, d2=dist2, d3=dist3, MoreArgs=list(theta=NULL, rad=1))
angles
}
