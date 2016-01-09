print.sphwin <- function(x, ...) {
	stopifnot(inherits(x, "sphwin"))
	cat("Window: ")
	switch(x$type,
		sphere = {
			cat("sphere with radius ", x$rad, fill=TRUE)
		},
		band = {
			if(x$param[1]==0) {
				cat("cap with centre (", x$ref, ") and ",
                                    "base congruent to colatitude ",
                                    paste0(x$param[2], "."),
                                    fill=TRUE)
			}
			else if(x$param[2]==pi && x$ref*(pi-x$ref) !=0) {
				cat("cap with centre ",
                                    c(pi-x$ref[1], (x$ref[2] + pi) %% (2*pi)),
                                    " and base congruent to ",
                                    "colatitude ",
                                    paste0(pi-x$param[1], "."),
                                    fill=TRUE)
			}
			else if(x$param[2]==pi && x$ref*(pi-x$ref) ==0) {
				cat("cap with centre ",
                                    c(pi-x$ref[1], 0),
                                    " and base congruent to ",
                                    "colatitude ",
                                    paste0(pi-x$param[1], "."),
                                    fill=TRUE)
			}
			else {
				cat("band with normal (", x$ref, ")",
                                    " and bounded by circles ",
                                    "congruent to colatitudes ",
                                    x$param[1], " and ", x$param[2], ".",
                                    fill=TRUE)
			}
		},
		bandcomp = {
			cat("the union of the cap centred at (", x$ref, ")",
                            " and base congruent to colatitude ",
                            paste0(x$param[1], ","),
                            " and the cap centred ",
                            "diametrically opposite",
                            " with base congruent to", " colatitude ",
                            paste0(x$param[2], "."),
                            fill=TRUE)
		},
		wedge = {
			cat("wedge with vertices (", x$ref, ") and (",
                            c(pi-x$ref[1], (x$ref[2]+pi)%%(2*pi)), "), ",
                            "dihedral angle ", x$param[1],
                            " and angular correction ",
                            paste0(x$param[2], "."),
                            fill=TRUE)
		},
		polygon = {
			cat("polygon with ", nrow(x$param), "vertices.",
                            fill=TRUE)
		},
		quadrangle = {
			cat("quadrangle formed by the",
                            " intersection of the band bounded by ",
                            "circles congruent to colatitude ", x$param[1],
                            " and ", x$param[2], ", and the wedge with",
                            " a vertex at (", x$ref, "), dihedral angle ",
                            x$param[3], " and angular correction ",
                            x$param[4], ".",
                            fill=TRUE)
		},
		{
			stop("Unsupported shape type")
		}
	)
        if(x$type != "sphere") 
          cat("Subset of sphere with radius ", x$rad, fill=TRUE)
}
