### =========================================================================
### GenomicRangesList objects
### -------------------------------------------------------------------------
###
### A List of GenomicRanges objects. Subclasses not necessarily have the same
### "compound" semantics as GRangesList.
###


setClass("GenomicRangesList",
    contains="List",
    representation(
        "VIRTUAL",
        elementMetadata="DataFrame"
    ),
    prototype(
        elementType="GenomicRanges"
    )
)

setClass("SimpleGenomicRangesList",
    contains=c("GenomicRangesList", "SimpleList")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GenomicRangesList <- function(...) {
  args <- list(...)
  if (length(args) == 1 && is.list(args[[1]]))
    args <- args[[1]]
  S4Vectors:::new_SimpleList_from_list("SimpleGenomicRangesList", args)
}
