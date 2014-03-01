### =========================================================================
### Generic functions for squeezing the ranges out of a range-based object
### -------------------------------------------------------------------------

### Extract the ranges as a GRanges object.
setGeneric("granges", signature="x",
    function(x, use.mcols=FALSE, ...) standardGeneric("granges")
)

### Extract the ranges as a GRangesList object.
setGeneric("grglist", signature="x",
    function(x, use.mcols=FALSE, ...) standardGeneric("grglist")
)

### Extract the ranges as a RangesList object.
### TODO: This one should probably be in IRanges together with ranges(), which
### is another range-squeezer.
### TODO: For consistency the ranges() generic should also get the 'use.mcols'
### arg with default to FALSE.
setGeneric("rglist", signature="x",
    function(x, use.mcols=FALSE, ...) standardGeneric("rglist")
)

