### =========================================================================
### RangesMapping methods
### -------------------------------------------------------------------------
###

setMethod("granges",  "RangesMapping", function(x) {
  GRanges(space(x), ranges(x))
})

setAs("RangesMapping", "GenomicRanges", function(from) {
  gr <- granges(from)
  values(gr) <- DataFrame(as.matrix(matching(from)))
  gr
})
