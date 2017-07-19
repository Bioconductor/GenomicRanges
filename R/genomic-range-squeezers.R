### =========================================================================
### Generic functions for squeezing the genomic ranges out of a range-based
### object
### -------------------------------------------------------------------------


### Extract the genomic ranges as a GRanges object.
setGeneric("granges", signature="x",
    function(x, use.names=TRUE, use.mcols=FALSE, ...)
        standardGeneric("granges")
)

### Extract the genomic ranges as a GRangesList object.
setGeneric("grglist", signature="x",
    function(x, use.names=TRUE, use.mcols=FALSE, ...)
        standardGeneric("grglist")
)

### Pairs method.
setMethod("grglist", "Pairs", function(x, use.names=TRUE, use.mcols=FALSE) {
              stopifnot(isTRUEorFALSE(use.mcols))
              grl <- zipup(granges(first(x)), granges(second(x)))
              if (!use.mcols) {
                  mcols(grl) <- NULL
              }
              grl
          })

