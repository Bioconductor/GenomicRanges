### =========================================================================
### Set operations
### -------------------------------------------------------------------------

setMethod("union", c("GRanges", "GRanges"),
    function(x, y)
    {
        elementMetadata(x) <- NULL
        elementMetadata(y) <- NULL
        reduce(c(x, y), drop.empty.ranges=TRUE)
    }
)

setMethod("intersect", c("GRanges", "GRanges"),
    function(x, y)
    {
        elementMetadata(x) <- NULL
        elementMetadata(y) <- NULL
        seqnames <- unique(c(levels(seqnames(x)), levels(seqnames(y))))
        seqlengths <- c(seqlengths(x), seqlengths(y))[seqnames]
        if (!identical(seqlengths(x), seqlengths))
            seqlengths(x) <- seqlengths
        if (!identical(seqlengths(y), seqlengths))
            seqlengths(y) <- seqlengths
        if (IRanges:::anyMissing(seqlengths)) {
            maxs <-
              sapply(IRanges:::newCompressedList("CompressedIntegerList",
                               unlistData = c(end(x), end(y)),
                               splitFactor = c(seqnames(x), seqnames(y))),
                     function(x) if (length(x) > 0) max(x) else NA_integer_)
            seqlengths[names(maxs)] <- maxs
        }
        setdiff(x, gaps(y, end = seqlengths))
    }
)

setMethod("setdiff", c("GRanges", "GRanges"),
    function(x, y) {
        elementMetadata(x) <- NULL
        elementMetadata(y) <- NULL
        seqnames <- unique(c(levels(seqnames(x)), levels(seqnames(y))))
        seqlengths <- c(seqlengths(x), seqlengths(y))[seqnames]
        if (!identical(seqlengths(x), seqlengths))
            seqlengths(x) <- seqlengths
        if (!identical(seqlengths(y), seqlengths))
            seqlengths(y) <- seqlengths
        if (IRanges:::anyMissing(seqlengths)) {
            maxs <-
              sapply(IRanges:::newCompressedList("CompressedIntegerList",
                               unlistData = c(end(x), end(y)),
                               splitFactor = c(seqnames(x), seqnames(y))),
                     function(x) if (length(x) > 0) max(x) else NA_integer_)
            seqlengths[names(maxs)] <- maxs
        }
        gaps(union(gaps(x, end = seqlengths), y), end = seqlengths)
    }
)


### =========================================================================
### Parallel set operations
### -------------------------------------------------------------------------

setMethod("punion", c("GRanges", "GRanges"),
    function(x, y, fill.gap = FALSE, ...)
    {
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        if (!setequal(levels(seqnames(x)), levels(seqnames(y))) ||
            !all((seqnames(x) == seqnames(y)) &
                 (strand(x) == strand(y))))
            stop("'x' and 'y' must elements have compatable 'seqnames' ",
                 "and 'strand' values")
        GRanges(seqnames(x),
                callGeneric(ranges(x), ranges(y), fill.gap = fill.gap),
                strand(x),
                seqlengths = seqlengths(x))
    }
)

setMethod("punion", c("GRangesList", "GRanges"),
    function(x, y, fill.gap = FALSE, ...)
    {
        n <- length(x)
        if (n != length(y)) 
            stop("'x' and 'y' must have the same length")
        elementMetadata(x@unlistData) <- NULL
        elementMetadata(y) <- NULL
        ans <-
          split(c(x@unlistData, y), 
                c(Rle(seq_len(n), elementLengths(x)), Rle(seq_len(n))))
        names(ans) <- names(x)
        ans
    }
)

setMethod("punion", c("GRanges", "GRangesList"),
    function(x, y, fill.gap = FALSE, ...)
    {
        callGeneric(y, x)
    }
)


setMethod("pintersect", c("GRanges", "GRanges"),
    function(x, y, resolve.empty = c("none", "max.start", "start.x"), ...)
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
        GRanges(seqnames(x), ansRanges, ansStrand, seqlengths = seqlengths(x))
    }
)

setMethod("pintersect", c("GRangesList", "GRanges"),
    function(x, y, resolve.empty = c("none", "max.start", "start.x"), ...)
    {
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
        if (ncol(elementMetadata(x@unlistData)) > 0)
            elementMetadata(x@unlistData) <- NULL
        if (ncol(elementMetadata(y)) > 0)
            elementMetadata(y) <- NULL
        x <- x[ok]
        y <- rep(y, sum(ok))
        x@unlistData@ranges <-
          callGeneric(x@unlistData@ranges, y@ranges,
                      resolve.empty = "start.x")
        x[width(x) > 0L]
    }
)

setMethod("pintersect", c("GRanges", "GRangesList"),
    function(x, y, resolve.empty = c("none", "max.start", "start.x"), ...)
    {
        callGeneric(y, x)
    }
)


setMethod("psetdiff", c("GRanges", "GRanges"),
    function(x, y, ...)
    {
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        if (!setequal(levels(seqnames(x)), levels(seqnames(y))))
            stop("'x' and 'y' must elements have compatable 'seqnames'")
        ok <-
          (seqnames(x) == seqnames(y)) & compatableStrand(strand(x), strand(y))
        ansRanges <- ranges(x)
        seqselect(ansRanges, ok) <-
          callGeneric(seqselect(ranges(x), ok), seqselect(ranges(y), ok))
        ansStrand <- strand(x)
        resolveStrand <- as(ansStrand == "*", "IRanges")
        if (length(resolveStrand) > 0)
            ansStrand[as.integer(resolveStrand)] <-
              seqselect(strand(y), resolveStrand)
        GRanges(seqnames(x), ansRanges, ansStrand, seqlengths = seqlengths(x))
    }
)

setMethod("psetdiff", c("GRanges", "GRangesList"),
    function(x, y, ...)
    {
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        if (!setequal(levels(seqnames(x)), levels(seqnames(y@unlistData))))
            stop("'x' and 'y' must elements have compatable 'seqnames'")
        ok <-
          (rep(seqnames(x), elementLengths(y)) == seqnames(y@unlistData)) &
           compatableStrand(rep(strand(x), elementLengths(y)), strand(y@unlistData))
        if (!all(ok)) {
            ok <-
              new2("CompressedLogicalList", unlistData = as.vector(ok),
                   partitioning = y@partitioning)
            y <- y[ok]
        }
        ansRanges <- gaps(ranges(y), start = start(x), end = end(x))
        ansSeqnames <- rep(seqnames(x), elementLengths(ansRanges))
        ansStrand <- rep(strand(x), elementLengths(ansRanges))
        ansGRanges <-
          GRanges(ansSeqnames, unlist(ansRanges, use.names = FALSE), ansStrand,
                  seqlengths = seqlengths(x))
        new2("GRangesList", elementMetadata = new("DataFrame", nrows = length(x)),
             unlistData = ansGRanges, partitioning = ansRanges@partitioning,
             check=FALSE)
    }
)
