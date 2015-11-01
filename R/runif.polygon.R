runif.polygon <- function(n, win) {
  stopifnot(inherits(win, "sphwin") && win$type=="polygon")
  esphere <- sphwin(type="sphere", rad=win$rad)
  output <- matrix(, nrow=0, ncol=2)
  while((need <- (n - nrow(output))) > 0) {
    proposed <- runif.sphere(need, win=esphere)
    accept <- in.W(points=proposed, win=win)
    if(any(accept)) 
      output <- rbind(output, proposed[accept, , drop=FALSE])
  }
  if(nrow(output) > n) 
    output <- output[seq_len(nrow), , drop=FALSE]
  return(output)
}
