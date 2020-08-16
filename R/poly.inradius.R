poly.inradius <- function(win, ncolat=100, nlon=100) {
	stopifnot(inherits(win,"sphwin") && win$type=="polygon")
	param <- win$param
	rad <- win$rad
	r <- seq(0,pi*rad,length=1024)
	gridsph <- gridmat.nlon(colats=range(param[,1]),
                                     lons=range(param[,2]), ncolat=100, nlon=100)
        ea <- polysph.area.Wr.grid(win=win, points=gridsph$gridrefs,
                                        nlon=gridsph$nlon, r=r)
	rref <- sum(ea > 0)
	inradius <- r[rref]
	inradius
}
        