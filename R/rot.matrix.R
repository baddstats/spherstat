rot.matrix <-
function(northpole, rad=1) {
if(class(northpole) != "matrix") {northpole1 <- t(matrix(northpole))} else {northpole1 <- northpole}
A <- matrix(nrow=3, ncol=3)
A[3,] <- convert3(northpole1)/rad
theta <- northpole1[1,1]
phi <- northpole1[1,2]
A[1,] <- c(cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta))
A[2,] <- c(-sin(phi), cos(phi), 0)
for(i in 1:3){for(j in 1:3) {if(-(10^(-16)) <= A[i,j] && (10^(-16)) >= A[i,j]) {A[i,j] <- 0}}}
A
}
