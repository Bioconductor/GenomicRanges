### =========================================================================
### DelegatingGenomicRanges objects
### -------------------------------------------------------------------------
###
### Virtual class that delegates GenomicRanges data access to a
### GenomicRanges delegate.
###

setClass("DelegatingGenomicRanges",
         representation(delegate="GenomicRanges"),
         contains=c("GenomicRanges", "VIRTUAL"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Slot getters and setters.
###

setMethod("seqnames", "DelegatingGenomicRanges",
          function(x) seqnames(x@delegate))

setMethod("ranges", "DelegatingGenomicRanges",
          function(x, ...) ranges(x@delegate, ...))

setMethod("strand", "DelegatingGenomicRanges", function(x) strand(x@delegate))
setMethod("seqinfo", "DelegatingGenomicRanges", function(x) seqinfo(x@delegate))

setMethod("update", "DelegatingGenomicRanges", function (object, ...) {
  object@delegate <- update(object@delegate, ...)
  object
})
