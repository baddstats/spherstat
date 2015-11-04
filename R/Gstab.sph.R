Gstab.sph <- function (X, win=sphwin(type="sphere"), r=NULL, ...) {
    G <- Gsphere(X=X, win=win, r=r, ...)
    rad <- X$win$rad
    Gstab <- eval.fv(asin(sqrt(G))/(pi/2))
    if (any(varcols <- colnames(G) %in% c("rip", "ls"))) {
        r <- with(Gstab, .x)
        Gstab[, varcols] <- as.data.frame(G)[, varcols]/(2 * pi * (rad^2)* (1-cos(r/rad)))
        n <- npoints(X)
        A <- area.sphwin(X$win)
        if (any(colnames(G) == "rip")) 
            Gstab[r == 0, "rip"] <- (2 * A/(n - 1)^2)/(4 * pi)
        if (any(colnames(G) == "ls")) 
            Gstab[r == 0, "ls"] <- (2 * A/(n * (n - 1)))/(4 * pi)
    }
    Gstab <- rebadge.fv(Gstab, quote(Gstab(r)), "Gstab", names(G), new.labl = attr(G, 
        "labl"))
    return(Gstab)
}
