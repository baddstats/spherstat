eroded.areas.sphwin <-
function(win=sphwin(type="sphere"), r=NULL, method="exact", ...) {
stopifnot(inherits(win, c("sp2", "sp3", "sphwin")))
if(inherits(win, "sphwin")) {win <- win} else {win <- win$win}
rad <- win$rad
if(is.null(r)) { r <- switch(win$type,
sphere = {seq(0, pi*rad, length=512)},
band = {seq(0, rad*sum(win$param)/2, length=512)},
bandcomp = {seq(0, rad*max(win$param[1], pi-win$param[2]), length=512)},
wedge = {ifelse(win$param[1] < pi, seq(0, rad*win$param[1], length=512), seq(0, rad*pi/2, length=512))},
quadrangle = {stop("r needs to be specified for this window")},
polygon = {stop("r needs to be specified for this window")},
{stop("Unrecognised window type")}
)}
eroded.area <- switch(win$type,
sphere = {ea <- rep(area.sphwin(w=win), times=512)
ea},
band = {ea <- c()
if(win$param[1]==0) {
if(max(r) > win$param[2]) {
stop("r range too wide for window")
} else {
ea <- 2*pi*(1-cos(win$param[2]-r/rad))
}} else {
if(win$param[2]==pi) {
if(max(r) > pi-win$param[1]) {
stop("r range too wide for window")
} else {
ea <- 2*pi*(1-cos((pi-win$param[1])-r/rad))
}} else {
if(max(r) > sround((win$param[2]-win$param[1])/2)) {
stop("r range too wide for window")
} else {
sbp1 <- win$param[1]+(r/rad)
sbp2 <- win$param[2]-(r/rad)
ea <- 2*pi*(rad^2)*(cos(sbp1)-cos(sbp2))
}}}
ea
},
bandcomp = {
if(sround(max(r) - max(win$param[1], pi-win$param[2]))>0) {
stop("r range too wide for window")
} else {
ea <- pmax(0, 2*pi*(1-cos(win$param[1]-r/rad))) + pmax(0, 2*pi*(1-cos((pi-win$param[2])-r/rad)))}
ea
},
wedge = {switch(method,
exact = {polysph.area.Wr.exact(win=win, r=r)},
integral = {polysph.area.Wr.int(win=win, r=r)},
grid = {gridsph <- gridmat.nlon(colats=c(0, pi), lons=c(0, win$param[1]), ...)
ea <- polysph.area.Wr.grid(win=win, points=gridsph$gridrefs, nlon=gridsph$nlon, r=r)
ea},
stop("method not recognised"))},
polygon = {
gridsph <- gridmat.nlon(colats=range(win$param[,1]), lons=range(win$param[,2]), ...)
ea <- polysph.area.Wr.grid(win=win, points=gridsph$gridrefs, nlon=gridsph$nlon, r=r)
ea},
quadrangle = {
gridsph <- gridmat.nlon(colats=(win$param[1:2]), lons=c(0, win$param[3]), ...)
ea <- polysph.area.Wr.grid(win=win, points=gridsph$gridrefs, nlon=gridsph$nlon, r=r)
ea},
{stop("Unrecognised window type")})
eroded.area <- as.numeric(eroded.area)
eroded.area
}
## if((win$param[1]==0 & max(r) > win$param[2]) || (win$param[2]==pi & max(r) > pi-win$param[1]) || (win$param!=0
## if(min(sround(((win$param[2]-win$param[1])/2)-r)) < 0) {