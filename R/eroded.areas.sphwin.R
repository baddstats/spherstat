eroded.areas.sphwin <- function(win=sphwin(type="sphere"),
                                r=NULL, method="exact", ...) {
  stopifnot(inherits(win, c("sp2", "sp3", "sphwin")))
  if(inherits(win, "sphwin")) {
    win <- win
  } else {
    win <- win$win
  }
  rad <- win$rad
  param <- win$param
  if(is.null(r)) {
    rmax <- rmax.rule.sphwin(win)
    r <- seq(0, rmax, length=512)
  }
  eroded.area <-
    switch(win$type,
           sphere = rep(area.sphwin(w=win), times=length(r)),
           band = {
             if(param[1]==0) {
               2*pi*(rad^2)*(1-cos(pmax(0, param[2]-r/rad)))
             } else if(param[2]==pi) {
               2*pi*(rad^2)*(1-cos(pmax(0, (pi-param[1])-r/rad)))
             } else {
               sbp1 <- param[1]+(r/rad)
               sbp2 <- param[2]-(r/rad)
               2*pi*(rad^2)* pmax(cos(sbp1)-cos(sbp2), 0)
             }
           },
           bandcomp = {
             2*pi*(rad^2)*(
               (1-cos(pmax(0, param[1]-r/rad))) +
               (1-cos(pmax(0, (pi-param[2])-r/rad)))
               )
           },
           wedge = {
             switch(method,
                    exact = {polysph.area.Wr.exact(win=win, r=r)},
                    integral = {polysph.area.Wr.int(win=win, r=r)},
                    grid = {
                      gridsph <- gridmat.nlon(colats=c(0, pi),
                                              lons=c(0, param[1]), ...)
                      ea <- polysph.area.Wr.grid(win=win,
                                                 points=gridsph$gridrefs,
                                                 nlon=gridsph$nlon, r=r)
                      ea
                    },
                    stop("method not recognised"))
           },
           polygon = {
             gridsph <- gridmat.nlon(colats=range(param[,1]),
                                     lons=range(param[,2]), ...)
             ea <- polysph.area.Wr.grid(win=win, points=gridsph$gridrefs,
                                        nlon=gridsph$nlon, r=r)
             ea
           },
           quadrangle = {
             gridsph <- gridmat.nlon(colats=(param[1:2]),
                                     lons=c(0, param[3]), ...)
             ea <- polysph.area.Wr.grid(win=win, points=gridsph$gridrefs,
                                        nlon=gridsph$nlon, r=r)
             ea
           },
           {stop("Unrecognised window type")}
           )
  return(as.numeric(eroded.area))
}


rmax.rule.sphwin <- function(win) {
  stopifnot(inherits(win, "sphwin"))
  param <- win$param
  r1max <- switch(win$type,
                   sphere = pi,
                   band = mean(param),
                   bandcomp = max(param[1], pi - param[2]),
                   wedge = if(param[1] < pi) param[1] else pi/2,
                   quadrangle = ,
                   polygon = stop("r needs to be specified for this window"),
                   stop("Unrecognised window type")
                   )
  return(r1max * win$rad)
}
