### =========================================================================
### RangesMapping methods
### -------------------------------------------------------------------------
###

setMethod("granges",  "RangesMapping", function(x, use.mcols=FALSE) {
  if (!identical(use.mcols, FALSE))
    stop("\"granges\" method for RangesMapping objects ",
         "does not support the 'use.mcols' argument")
  GRanges(space(x), ranges(x))
})

setAs("RangesMapping", "GenomicRanges", function(from) {
  gr <- granges(from)
  mcols(gr) <- DataFrame(as.matrix(hits(from)))
  gr
})
