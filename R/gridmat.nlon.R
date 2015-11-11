gridmat.nlon <- function(colats, lons, ncolat=100, nlon=100) {
	stopifnot(colats[1] < colats[2] && lons[1] < lons[2])
	ctheta <- seq(from=cos(colats[1]), to=cos(colats[2]), length=ncolat+1)
	phi <- seq(from=lons[1], to=lons[2], length=nlon+1)
	thetas <- acos(ctheta[1:ncolat] + (ctheta[2:(ncolat+1)] - ctheta[1:ncolat])/2)
	phis <- phi[1:nlon]+((phi[2:(nlon+1)] - phi[1:nlon])/2)
	thetaslist <- rep(thetas, each=nlon)
	phislist <- rep(phis, ncolat)
	gridrefs <- matrix(cbind(thetaslist, phislist), ncol=2)
	gridmat <- list(gridrefs=gridrefs, nlon=nlon)
	gridmat
}
