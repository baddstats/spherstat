Fstab.sph <- function (X, refpoints, win=sphwin(type="sphere"), r=NULL, ...) {
    FF <- Fsphere(X, refpoints=refpoints, win=win, r=r, ...)
    rad <- X$win$rad
    Fstab <- eval.fv(asin(sqrt(FF))/(pi/2))
    if (any(varcols <- colnames(FF) %in% c("rip", "ls"))) {
        r <- with(Fstab, .x)
        Fstab[, varcols] <- as.data.frame(FF)[, varcols]/(2 * pi * (rad^2)* (1-cos(r/rad)))
        n <- npoints(X)
        A <- area.sphwin(X$win)
        if (any(colnames(FF) == "rip")) 
            Fstab[r == 0, "rip"] <- (2 * A/(n - 1)^2)/(4 * pi)
        if (any(colnames(FF) == "ls")) 
            Fstab[r == 0, "ls"] <- (2 * A/(n * (n - 1)))/(4 * pi)
    }
    Fstab <- rebadge.fv(Fstab, quote(Fstab(r)), "Fstab", names(FF), new.labl = attr(FF, 
        "labl"))
    return(Fstab)
}
