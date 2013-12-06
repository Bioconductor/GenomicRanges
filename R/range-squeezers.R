### =========================================================================
### Generic functions for squeezing the ranges out of a range-based object
### -------------------------------------------------------------------------

### Extract the ranges as a GRanges object.
setGeneric("granges", function(x, ...) standardGeneric("granges"))

### Extract the ranges as a GRangesList object.
setGeneric("grglist", function(x, ...) standardGeneric("grglist"))

### Extract the ranges as a RangesList object.
### TODO: Maybe this one should be in IRanges?
setGeneric("rglist", function(x, ...) standardGeneric("rglist"))

