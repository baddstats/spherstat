cround <- function(x) {
	if(inherits(x, "list") || inherits(x, "data.frame")) {
		for(i in 1:length(x)) {
			x[[i]] <- ifelse(x[[i]] < -1, -1, ifelse(x[[i]] > 1, 1, x))
		}
	}
	else {
		ifelse(x < -1, -1, ifelse(x > 1, 1, x))
	}
}