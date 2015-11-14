print.sphwin <- function(x, ...) {
	stopifnot(inherits(x, "sphwin"))
	cat("Window: ")
	switch(x$type,
		sphere = {
			cat("sphere with radius ", x$rad, "\n")
		},
		band = {
			if(x$param[1]==0) {
				cat("cap with centre (", x$ref, ") and base congruent to colatitude ", x$param[2], ".  Subset of sphere with radius ", x$rad, ".", "\n")
			}
			else if(x$param[2]==pi && x$ref*(pi-x$ref) !=0) {
				cat("cap with centre ", c(pi-x$ref[1], (x$ref[2] + pi) %% (2*pi)), " and base congruent to colatitude ", pi-x$param[1], ".  Subset of sphere with radius ", x$rad, ".", "\n")
			}
			else if(x$param[2]==pi && x$ref*(pi-x$ref) ==0) {
				cat("cap with centre ", c(pi-x$ref[1], 0), " and base congruent to colatitude ", pi-x$param[1], ".  Subset of sphere with radius ", x$rad, ".", "\n")
			}
			else {
				cat("band with normal (", x$ref, ") and bounded by circles congruent to colatitudes ", x$param[1], " and ", x$param[2], ". Subset of sphere with radius ", x$rad, ".", "\n")
			}
		},
		bandcomp = {
			cat("the union of the cap centred at (", x$ref, ") and base congruent to colatitude ", x$param[1], ", and the cap centred diametrically opposite with base congruent to colatitude ", x$param[2], ".  Subset of sphere with radius ", x$rad, ".", "\n")
		},
		wedge = {
			cat("wedge with vertices (", x$ref, ") and (", c(pi-x$ref[1], (x$ref[2]+pi)%%(2*pi)), "), dihedral angle ", x$param[1], " and angular correction ", x$param[2], ".  Subset of sphere with radius ", x$rad, ".", "\n")
		},
		polygon = {
			cat("polygon with ", nrow(x$param), "vertices.  Subset of sphere with radius ", x$rad, ".", "\n")
		},
		quadrangle = {
			cat("quadrangle formed by the intersection of the band bounded by circles congruent to colatitude ", x$param[1], " and ", x$param[2], ", and the wedge with a vertex at (", x$ref, "), dihedral angle ", x$param[3], " and angular correction ", x$param[4], ".  Subset of sphere with radius ", x$rad, ".", "\n")
		},
		{
			stop("Unsupported shape type")
		}
	)
}