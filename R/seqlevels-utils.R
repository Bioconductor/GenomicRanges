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

setMethod("keepSeqlevels",  c("GenomicRanges", "GappedAlignments"),
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

setMethod("keepSeqlevels",  c("GRangesList", "GappedAlignments"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GenomicRanges", "character"),
            function(x, value, ...)
{
    if (!any(value %in% seqlevels(x)))
        stop("none of the values in 'value' are present in 'x'")
    suppressWarnings(seqlevels(x, force=TRUE) <- value) 
    x
})

setMethod("keepSeqlevels",  c("GRangesList", "character"),
            function(x, value, ...)
{
    if (!any(value %in% seqlevels(x)))
        stop("none of the values are present in 'x'")

    grlSeqnames <- seqnames(unlist(x))
    idx <- rep(runValue(grlSeqnames) %in% value, runLength(grlSeqnames))
    grReduced <- GenomicRanges:::deconstructGRLintoGR(x)[idx]
    grlReduced <- GenomicRanges:::reconstructGRLfromGR(grReduced, x)

    seqlevels(grlReduced) <- seqlevels(x)[seqlevels(x) %in% value]
    metadata(grlReduced) <- metadata(x)
    grlReduced[elementLengths(grlReduced) != 0] 
})

setMethod("keepSeqlevels",  c("GappedAlignments", "GenomicRanges"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GappedAlignments", "GRangesList"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GappedAlignments", "GappedAlignments"),
            function(x, value, ...)
{
    value <- seqlevels(value)
    callGeneric(x, value, ...)
})

setMethod("keepSeqlevels",  c("GappedAlignments", "character"),
            function(x, value, ...)
{
    if (!any(value %in% seqlevels(x)))
        stop("none of the values in 'value' are present in 'x'")
    suppressWarnings(seqlevels(x, force=TRUE) <- value) 
    x
})


## renameSeqlevels : 

.renameSeqlevels <- function(x, value, ...)
{
    if (is.null(names(value)))
        stop("elements in 'value' must be named")
    old <- names(value)
    new <- unlist(value, use.names=FALSE)
    if (!any(old %in% seqlevels(x)))
        stop("none of the values in 'value' are present in 'x'")
    seqlevels(x)[seqlevels(x) == old] <- new
    x 
}

setGeneric("renameSeqlevels", signature = c("x", "value"),
    function(x, value, ...)
    standardGeneric("renameSeqlevels")
)

setMethod("renameSeqlevels",  c(x="GappedAlignments", value="character"),
    function(x, value, ...) .renameSeqlevels(x, value, ...)
)

setMethod("renameSeqlevels",  c("GRangesList", "character"),
    function(x, value, ...) .renameSeqlevels(x, value, ...)
)

setMethod("renameSeqlevels",  c("GenomicRanges", "character"),
    function(x, value, ...) .renameSeqlevels(x, value, ...)
)

