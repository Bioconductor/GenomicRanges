### =========================================================================
### RangesMapping methods
### -------------------------------------------------------------------------
###

setMethod("granges",  "RangesMapping", function(x) {
  GRanges(space(x), ranges(x))
})

setAs("RangesMapping", "GenomicRanges", function(from) {
  gr <- granges(from)
  mcols(gr) <- DataFrame(as.matrix(hits(from)))
  gr
})
