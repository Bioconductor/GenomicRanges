### =========================================================================
### GenomicRangesList objects
### -------------------------------------------------------------------------
###
### A SimpleList of GenomicRanges objects. Does NOT have the same
### "compound" semantics as GRangesList.
###

setClass("GenomicRangesList",
         prototype = prototype(elementType = "GenomicRanges"),
         contains = "SimpleList")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GenomicRangesList <- function(...) {
  args <- list(...)
  if (length(args) == 1 && is.list(args[[1]]))
    args <- args[[1]]
  IRanges:::newSimpleList("GenomicRangesList", args)
}

