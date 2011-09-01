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
    x <- x[seqnames(x) %in% value]
    seqlevels(x) <- seqlevels(x)[seqlevels(x) %in% value]
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
    x <- x[seqnames(x) %in% value]
    seqlevels(x) <- seqlevels(x)[seqlevels(x) %in% value]
    x
})



## renameSeqlevels : 

.renameSeqlevels <- function(x, value, ...)
{
    old <- names(value)
    new <- unlist(value, use.names=FALSE)
    if (!any(old %in% seqlevels(x)))
        stop("none of the values in 'value' are present in 'x'")
    seqlevels(x)[seqlevels(x) %in% old] <- new
    x 
}

setGeneric("renameSeqlevels", signature = c("x", "value"),
    function(x, value, ...)
    standardGeneric("renameSeqlevels")
)

setMethod("renameSeqlevels",  c(x="GappedAlignments", value="list"),
    function(x, value, ...) .renameSeqlevels(x, value, ...)
)

setMethod("renameSeqlevels",  c("GRangesList", "list"),
    function(x, value, ...) .renameSeqlevels(x, value, ...)
)

setMethod("renameSeqlevels",  c("GenomicRanges", "list"),
    function(x, value, ...) .renameSeqlevels(x, value, ...)
)

.SubsetGRListAtRangesLevel <- function(grl, idx, ...)
{
    if (is.character(idx)) {
        orig <- idx
        idx <- match(idx, names(grl@unlistData))
    }
    grReduced <-
        GenomicRanges:::deconstructGRLintoGR(grl)[idx,,drop=FALSE]
    grlReduced <-
        GenomicRanges:::reconstructGRLfromGR(grReduced, grl)
    grlReduced[elementLengths(grlReduced) != 0]
}


