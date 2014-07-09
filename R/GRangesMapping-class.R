### =========================================================================
### GRangesMapping objects
### -------------------------------------------------------------------------
###
### A GRangesMapping encodes a mapping of a set of ranges to some other
### coordinate space.
###

## Conceptually, a GRangesMapping is a matching of each query range to
## one or more elements in a subject. The geometry of the query range
## is then transformed according to that matching. Thus, this data
## class combines a Hits object with a set of transformed ranges.

setClass("GRangesMapping",
         representation(hits = "Hits", granges = "GRanges"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

hits <- function(x) slot(x, "hits")

setMethod("granges",  "GRangesMapping", function(x, use.mcols=FALSE) {
  if (!identical(use.mcols, FALSE))
    stop("\"granges\" method for RangesMapping objects ",
         "does not support the 'use.mcols' argument")
  slot(x, "granges")
})

setMethod("dim", "GRangesMapping", function(x) dim(hits(x)))

setMethod("length", "GRangesMapping", function(x) length(granges(x)))

setMethod("subjectHits", "GRangesMapping", function(x) subjectHits(hits(x)))

setMethod("queryHits", "GRangesMapping", function(x) queryHits(hits(x)))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

setAs("GRangesMapping", "RangedData", function(from) {
  RangedData(ranges(granges(from)), space = seqnames(granges(from)), 
             as(hits(from), "DataFrame"))
})


setAs("GRangesMapping", "GenomicRanges", function(from) {
  gr <- granges(from)
  mcols(gr) <- DataFrame(as.matrix(hits(from)))
  gr
})
