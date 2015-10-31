runif.polygon <-
function(n, win) {
stopifnot(inherits(win, "sphwin") && win$type=="polygon")
rad <- win$rad
# ajb: 'rad' is defined but not used
# ajb: Query correctness of code when rad != 1
output <- t(matrix(0, nrow=2, ncol=1))
while(nrow(output) < (n+1)) {
point.test <- t(matrix(runif.sphere(1, win=sphwin(type="sphere", rad=1))))
if(in.W(points=point.test, win=win)) {
output <- rbind(output, point.test)
}
}
output <- output[2:n+1,]
output
}
