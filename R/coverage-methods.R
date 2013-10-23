### =========================================================================
### coverage methods
### -------------------------------------------------------------------------
###

setMethod("coverage", "GenomicRanges",
    function(x, shift=0L, width=NULL, weight=1L,
                method=c("auto", "sort", "hash"))
    {
        if (isSingleString(weight)) {
            x_mcols <- mcols(x)
            if (!is(x_mcols, "DataTable")
             || sum(colnames(x_mcols) == weight) != 1L)
                stop("'mcols(x)' has 0 or more than 1 \"",
                     weight, "\" columns")
            weight <- x_mcols[[weight]]
        }
        y <- split(ranges(x), seqnames(x))
        circle.length <- seqlengths(x)
        circle.length[!(isCircular(x) %in% TRUE)] <- NA_integer_
        IRanges:::.CompressedIRangesList.coverage(y,
                                        shift=shift,
                                        width=width,
                                        weight=weight,
                                        circle.length=circle.length,
                                        method=method)
    }
)

setMethod("coverage", "GRangesList",
    function(x, shift=0L, width=NULL, weight=1L, ...)
    {
        coverage(x@unlistData, shift=shift, width=width, weight=weight, ...)
    }
)

setMethod("coverage", "GAlignments",
    function(x, shift=0L, width=NULL, weight=1L, drop.D.ranges = FALSE, ...)
        coverage(grglist(x, drop.D.ranges = drop.D.ranges),
                 shift=shift, width=width, weight=weight, ...)
)

setMethod("coverage", "GAlignmentPairs",
    function(x, shift=0L, width=NULL, weight=1L, drop.D.ranges = FALSE, ...)
          coverage(grglist(x, drop.D.ranges = drop.D.ranges),
                   shift=shift, width=width, weight=weight, ...)
)

