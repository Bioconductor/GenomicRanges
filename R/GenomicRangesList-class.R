### =========================================================================
### GenomicRangesList objects
### -------------------------------------------------------------------------
###
### A List of GenomicRanges objects. Subclasses not necessarily have the same
### "compound" semantics as GRangesList.
###

setClass("GenomicRangesList",
         prototype = prototype(elementType = "GenomicRanges"),
         contains = "List")

setClass("SimpleGenomicRangesList",
         contains = c("GenomicRangesList", "SimpleList"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GenomicRangesList <- function(...) {
  args <- list(...)
  if (length(args) == 1 && is.list(args[[1]]))
    args <- args[[1]]
  IRanges:::newSimpleList("SimpleGenomicRangesList", args)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("RangedDataList", "GenomicRangesList",
      function(from) GenomicRangesList(lapply(from, as, "GRanges")))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utilities.
###

setMethod("stack", "GenomicRangesList", function(x, indName = "sample") {
  x_flat <- unlist(x, use.names = FALSE)
  values(x_flat) <- cbind(IRanges:::.stack.ind(x, indName), values(x_flat))
  x_flat
})
