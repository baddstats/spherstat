rMatClust.sp2 <- function(kappa, r, mu, win, parents=FALSE) {
stopifnot(inherits(win, "sphwin")&& kappa > 0 && mu > 0 && r > 0)
rad <- win$rad
rp <- rpoispp.sp2(lambda=kappa, win=win, as.sp=FALSE)
rpl <- nrow(rp)
rM <- rp
for(i in 1:rpl) {
nlam<-rpois(1, mu)
if(nlam > 0) {
daughtwin <- sphwin(type="band", param=c(0, r/rad), ref=rp[i,])
rMat1 <-  runif.sphwin(n=nlam,  win=daughtwin, as.sp=FALSE)
inrm <- in.W(points=rMat1, win=win)
for(j in 1:length(inrm)) {
if(inrm[j]==TRUE) {rM <- rbind(rM, rMat1[j,])
}}}
}
if(!parents) {output <- rM[(rpl+1):nrow(rM),]}  else {output <- rM}
out1 <- sp2(X=output, win=win)
out1
}
