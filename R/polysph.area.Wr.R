polysph.area.Wr <-
function(method="grid", ...) {
stopifnot(method=="integral" || method=="grid" || method=="exact")
areaWr <- switch(method,
  		"grid"={polysph.area.Wr.grid(...)},
		"integral"= {polysph.area.Wr.int(...)},
		"exact"={polysph.area.Wr.exact(...)}
		)
areaWr
}
