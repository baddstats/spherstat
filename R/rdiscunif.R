rdiscunif <- function(n, min, max, steps=1) {
m <- (max-min)/steps
stopifnot(floor(m)==m && min <= max)
output1 <- runif(n, 0, 1)
output <- (ceiling(output1*(m+1))-1)*steps+min
output
}
