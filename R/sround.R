sround <- function(x) {
	if(inherits(x, "list") || inherits(x, "data.frame")) {
		for(i in 1:length(x)) {
			ifelse(abs(x[[i]]) < 1e-7, 0, x[[i]])
		}
	}
	else {
			ifelse(abs(x) < 1e-7 , 0, x)
	}
}

