Jsphere <- function(X, refpoints=NULL, r=NULL, ..., correction = NULL) 
{
    W <- X$win
    FF <- Fsphere(X, refpoints=refpoints, win=W, r=r, ...)
    r <- FF$r
    G <- Gsphere(X, win=W, r=r, ...)
    rmax <- max(r)
    Z <- fv(data.frame(r = r, theo = 1), "r", quote(J(r)), 
            "theo", . ~ r, c(0, rmax), c("r", "%s[pois](r)"), 
            c("distance argument r", "theoretical Poisson %s"), fname = "J")
    ratio <- function(a, b) {
        result <- a/b
        result[b == 0] <- NA
        result
    }
    Fnames <- names(FF)
    Gnames <- names(G)
    if ("raw" %in% Gnames && "raw" %in% Fnames) {
        Jraw <- ratio(1 - G$raw, 1 - FF$raw)
        Z <- bind.fv(Z, data.frame(raw = Jraw), "hat(%s)[raw](r)", 
            "uncorrected estimate of %s", "raw")
        attr(Z, "alim") <- range(r[FF$raw <= 0.9])
    }
    if ("rs" %in% Gnames && "rs" %in% Fnames) {
        Jrs <- ratio(1 - G$rs, 1 - FF$rs)
        Z <- bind.fv(Z, data.frame(rs = Jrs), "hat(%s)[rs](r)", 
            "border corrected estimate of %s", "rs")
        attr(Z, "alim") <- range(r[FF$rs <= 0.9])
    }
    if ("han" %in% Gnames && "cs" %in% Fnames) {
        Jhan <- ratio(1 - G$han, 1 - FF$han)
        Z <- bind.fv(Z, data.frame(han = Jhan), "hat(%s)[han](r)", 
            "Hanisch-style estimate of %s", "han")
        attr(Z, "alim") <- range(r[FF$han <= 0.9])
    }
    if ("km" %in% Gnames && "km" %in% Fnames) {
        Jkm <- ratio(1 - G$km, 1 - FF$km)
        Z <- bind.fv(Z, data.frame(km = Jkm), "hat(%s)[km](r)", 
            "Kaplan-Meier estimate of %s", "km")
        attr(Z, "alim") <- range(r[FF$km <= 0.9])
    }
    if ("hazard" %in% Gnames && "hazard" %in% Fnames) {
        Jhaz <- G$hazard - FF$hazard
        Z <- bind.fv(Z, data.frame(hazard = Jhaz), "hazard(r)", 
            "Kaplan-Meier estimate of derivative of log(%s)")
    }
    nama <- names(Z)
    fvnames(Z, ".") <- rev(nama[!(nama %in% c("r", "hazard"))])
    attr(Z, "F") <- FF
    attr(Z, "G") <- G
    unitname(Z) <- unitname(X)
    return(Z)
}
