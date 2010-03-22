### =========================================================================
### pintersect methods
### -------------------------------------------------------------------------

setMethod("pintersect", c("GRanges", "GRanges"),
    function(x, y, resolve.empty=c("none", "max.start", "start.x"), ...)
    {
        resolve.empty <- match.arg(resolve.empty)
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        if (!setequal(levels(seqnames(x)), levels(seqnames(y))) ||
            !all((seqnames(x) == seqnames(y)) &
                  compatableStrand(strand(x), strand(y))))
            stop("'x' and 'y' must elements have compatable 'seqnames' ",
                 "and 'strand' values")
        ansRanges <-
          callGeneric(ranges(x), ranges(y), resolve.empty = resolve.empty)
        ansStrand <- strand(x)
        resolveStrand <- as(ansStrand == "*", "IRanges")
        if (length(resolveStrand) > 0)
            ansStrand[as.integer(resolveStrand)] <-
              seqselect(strand(y), resolveStrand)
        GRanges(seqnames(x), ansRanges, ansStrand)
    }
)

setMethod("pintersect", c("GRangesList", "GRanges"),
    function(x, y, resolve.empty=c("none", "max.start", "start.x"), ...)
    {
        resolve.empty <- match.arg(resolve.empty)
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        if (!setequal(levels(seqnames(x@unlistData)), levels(seqnames(y))))
            stop("'x' and 'y' must elements have compatable 'seqnames' ",
                 "and 'strand' values")
        ok <-
          (seqnames(x@unlistData) == rep(seqnames(y), elementLengths(x))) &
          compatableStrand(strand(x@unlistData), rep(strand(y), elementLengths(x)))
        ok <-
          new2("CompressedLogicalList", unlistData = as.vector(ok),
               partitioning = x@partitioning)
        if (ncol(x@unlistData) > 0)
            elementMetadata(x@unlistData) <- NULL
        if (ncol(y) > 0)
            elementMetadata(y) <- NULL
        x <- x[ok]
        y <- rep(y, sum(ok))
        x@unlistData@ranges <- callGeneric(x@unlistData@ranges, y@ranges)
        x[width(x) > 0L]
    }
)

setMethod("pintersect", c("GRanges", "GRangesList"),
    function(x, y, resolve.empty=c("none", "max.start", "start.x"), ...)
    {
        resolve.empty <- match.arg(resolve.empty)
        callGeneric(y, x, resolve.empty=resolve.empty)
    }
)
