Hstab.sph <- function (X, refpoints, win=sphwin(type="sphere"), rvals=seq(from=0, to=pi, length=512), ...) {
    HH <- Hsphere(X, refpoints=refpoints, win=win, rvals=rvals, ...)
    rad <- X$win$rad
    Hstab <- eval.fv(asin(sqrt(HH))/(pi/2))
    if (any(varcols <- colnames(HH) %in% c("rip", "ls"))) {
        r <- with(Hstab, .x)
        Hstab[, varcols] <- as.data.frame(HH)[, varcols]/(2 * pi * (rad^2)* (1-cos(r/rad)))
        n <- npoints(X)
        A <- area.sphwin(X$win)
        if (any(colnames(HH) == "rip")) 
            Hstab[r == 0, "rip"] <- (2 * A/(n - 1)^2)/(4 * pi)
        if (any(colnames(HH) == "ls")) 
            Hstab[r == 0, "ls"] <- (2 * A/(n * (n - 1)))/(4 * pi)
    }
    Hstab <- rebadge.fv(Hstab, quote(Hstab(r)), "Hstab", names(HH), new.labl = attr(HH, 
        "labl"))
    return(Hstab)
}
