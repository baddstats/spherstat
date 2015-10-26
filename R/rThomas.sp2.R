rThomas.sp2 <- function(kappa, sigma, mu, win, parents=FALSE) {
stopifnot(inherits(win, "sphwin") && kappa > 0 && sigma > 0 && mu > 0)
rp <- rpoispp.sp2(lambda=kappa, win=win, as.sp=FALSE)
rpl <- nrow(rp)
rThom2 <- rp
for(i in 1:rpl) {
np <- rpois(1, mu)
if(np > 0) {
rThom1 <-  rFisher(n=np, mode=rp[i,], kappa=sigma^-1, win=sphwin(type="sphere", rad=win$rad))$X
inrt <- in.W(points=rThom1, win=win)
for(j in 1:length(inrt)) {
if(inrt[j]) {rThom2 <- rbind(rThom2, rThom1[j,])
}}}}
if(!parents) {output <- rThom2[(rpl+1):nrow(rThom2),]}  else {output <- rThom2}
output <- sp2(X=output, win=win)
output
}
