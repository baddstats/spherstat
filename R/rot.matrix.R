rot.matrix <- function(northpole, rad=1) {
	if(!is.matrix(northpole)) {
		northpole <- t(matrix(northpole))
	}
	A <- matrix(nrow=3, ncol=3)
	theta <- northpole[1,1]
	phi <- northpole[1,2]
	ctheta <- cos(theta)
	stheta <- sin(theta)
	cphi <- cos(phi)
	sphi <- sin(phi)
	A[1,] <- c(ctheta*cphi, ctheta*sphi, -stheta)
	A[2,] <- c(-sphi, cphi, 0)
	A[3,] <- c(stheta*cphi, stheta*sphi, ctheta)
	A <- sround(cround(A))
	A
}
