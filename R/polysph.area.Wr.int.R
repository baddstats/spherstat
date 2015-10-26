polysph.area.Wr.int <-
function(win, r=0) {
stopifnot(inherits(win, c("sp2", "sp3", "sphwin")))
if(inherits(win, "sphwin")) {win <- win} else {win <- win$win}
stopifnot(win$type=="polygon" || win$type == "wedge")
rad <- win$rad
area.W <- area.sphwin(w=win)
lr <- length(r)
aW <- switch(win$type,
wedge = {
phi <- pi/2 - r/rad
theta <- win$param[1]/2
h <- pi - 2*sphsin(d1=r, d2=NULL, theta1=theta, theta2=pi/2, rad=rad)
psi0 <- 2*sphcos(d1=phi, d2=phi, d3=h, theta=NULL, rad=1)*rad
csp <- cumsum(rev(psi0))
area.Wr <- rev(c(0, csp[1:(lr-1)]*(r[2]-r[1])))*(rad^2)
area.Wr*(area.W/max(area.Wr))
},
stop("Unsupported window type"))
aW
}