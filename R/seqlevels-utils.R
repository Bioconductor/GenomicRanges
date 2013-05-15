### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Utility functions for subsetting and renaming seqlevels 
###

## keepSeqlevels :

setGeneric("keepSeqlevels", signature = c("x", "value"),
           function(x, value, ...)
           standardGeneric("keepSeqlevels")
)

setMethod("keepSeqlevels",  c("GenomicRanges", "GenomicRanges"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GenomicRanges", "GRangesList"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GenomicRanges", "GAlignments"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GRangesList", "GenomicRanges"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GRangesList", "GRangesList"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GRangesList", "GAlignments"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GenomicRanges", "character"),
            function(x, value, ...)
{
    .Deprecated("seqlevels")
    if (!identical(character(0), seqlevels(x)) &&
        !identical(character(0), value)) {
        str <- !value %in% seqlevels(x)
        if (any(str))
            warning("seqlevels ", paste(value[str], collapse=", "),
                    " not found in 'x'")
        if (all(str)) {
            x
        } else {
            suppressWarnings(seqlevels(x, force=TRUE) <- value[str == FALSE])
            x
        }
    } else {
        x
    }
})

setMethod("keepSeqlevels",  c("GRangesList", "character"),
            function(x, value, ...)
{
    .Deprecated("seqlevels")
    if (!identical(character(0), seqlevels(x)) &&
        !identical(character(0), value)) {
        str <- !value %in% seqlevels(x)
        if (any(str))
            warning("seqlevels ", paste(value[str], collapse=", "),
                    " not found in 'x'")
        if (all(str)) {
            x
        } else {
            grlSeqnames <- unlist(seqnames(x), use.names=FALSE)
            idx <- rep(runValue(grlSeqnames) %in% value, runLength(grlSeqnames))
            grReduced <- deconstructGRLintoGR(x)[idx]
            grlReduced <- reconstructGRLfromGR(grReduced, x)

            seqlevels(grlReduced) <- seqlevels(x)[seqlevels(x) %in% value]
            grlReduced[elementLengths(grlReduced) != 0] 
        }
    } else {
        x
    }
})

setMethod("keepSeqlevels",  c("GAlignments", "GenomicRanges"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GAlignments", "GRangesList"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GAlignments", "GAlignments"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GAlignments", "character"),
            function(x, value, ...)
{
    .Deprecated("seqlevels")
    if (!identical(character(0), seqlevels(x)) &&
        !identical(character(0), value)) {
        str <- !value %in% seqlevels(x)
        if (any(str))
            warning("seqlevels ", paste(value[str], collapse=", "), 
                    " not found in 'x'")
        if (all(str)) {
            x
        } else {
            suppressWarnings(seqlevels(x, force=TRUE) <- value[str == FALSE])
            x
        } 
    } else { 
        x
    }
})


## renameSeqlevels : 

.renameSeqlevels <- function(x, value, ...)
{
    .Deprecated("seqlevels")
    if (identical(character(0), seqlevels(x)))
        return(x)
    if (is.null(old <- names(value)))
        stop("'value' must be a named character vector")
    new <- unlist(value, use.names=FALSE)

    str <- !old %in% seqlevels(x)
    if (any(str)) {
        warning("seqlevels ", paste(old[str], collapse=", "), 
                " not found in 'x'")
    }
    ord <- match(old[str == FALSE], seqlevels(x))
    seqlevels(x)[ord] <- new[str == FALSE]
    x 
}

setGeneric("renameSeqlevels", signature = c("x", "value"),
    function(x, value, ...)
    standardGeneric("renameSeqlevels")
)

setMethod("renameSeqlevels",  c(x="GAlignments", value="character"),
    function(x, value, ...) .renameSeqlevels(x, value, ...)
)

setMethod("renameSeqlevels",  c("GRangesList", "character"),
    function(x, value, ...) .renameSeqlevels(x, value, ...)
)

setMethod("renameSeqlevels",  c("GenomicRanges", "character"),
    function(x, value, ...) .renameSeqlevels(x, value, ...)
)

