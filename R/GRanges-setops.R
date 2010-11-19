### =========================================================================
### Set operations
### -------------------------------------------------------------------------

### TODO: What's the impact of circularity on the set operations?

### NOTE: Doesn't properly merge 'seqinfo(x)' with 'seqinfo(y)' because
### 'c(x, y)' currently doesn't do it either. However, this should
### automatically get fixed when the "c" method itself gets fixed.
setMethod("union", c("GRanges", "GRanges"),
    function(x, y)
    {
        values(x) <- values(y) <- NULL  # so we can do 'c(x, y)' below
        reduce(c(x, y), drop.empty.ranges=TRUE)
    }
)

### 'x' must be a GRanges object.
### Returns a named integer vector where the names are guaranteed
### to be 'levels(seqnames(x))'.
.maxEndPerSequence <- function(x)
{
    ends_list <- IRanges:::newCompressedList("CompressedIntegerList",
                    unlistData=end(x),
                    splitFactor=seqnames(x))
    sapply(ends_list,
           function(ends)
               if (length(ends) > 0L) max(ends) else NA_integer_)
}

setMethod("intersect", c("GRanges", "GRanges"),
    function(x, y)
    {
        seqinfo(x) <- merge(seqinfo(x), seqinfo(y))
        seqinfo(y) <- suppressWarnings(merge(seqinfo(y), seqinfo(x)))
        seqlengths <- seqlengths(x)
        ## If the length of a sequence is unknown (NA), then we use
        ## the max end value found on that sequence in 'x' or 'y'.
        values(x) <- values(y) <- NULL  # so we can do 'c(x, y)' below
        seqlengths[is.na(seqlengths)] <-
            .maxEndPerSequence(c(x, y))[is.na(seqlengths)]
        setdiff(x, gaps(y, end = seqlengths))
    }
)

setMethod("setdiff", c("GRanges", "GRanges"),
    function(x, y)
    {
        seqinfo(x) <- merge(seqinfo(x), seqinfo(y))
        seqinfo(y) <- suppressWarnings(merge(seqinfo(y), seqinfo(x)))
        seqlengths <- seqlengths(x)
        ## If the length of a sequence is unknown (NA), then we use
        ## the max end value found on that sequence in 'x' or 'y'.
        values(x) <- values(y) <- NULL  # so we can do 'c(x, y)' below
        seqlengths[is.na(seqlengths)] <-
            .maxEndPerSequence(c(x, y))[is.na(seqlengths)]
        gaps(union(gaps(x, end = seqlengths), y), end = seqlengths)
    }
)


### =========================================================================
### Parallel set operations
### -------------------------------------------------------------------------

### FIXME: Why aren't all the parallel set operations using the same code
### for checking strand compatibility? E.g. "pintersect" and "psetdiff" use
### compatableStrand() for this but not "punion".

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
        ## TODO: 'seqinfo(x)' and 'seqinfo(y)' need to be merged. See TODO
        ## in "intersect" method at the top of this file for more details.
        ## Then the setequal check on the seqnames levels above won't be
        ## necessary anymore.
        GRanges(seqnames(x),
                callGeneric(ranges(x), ranges(y), fill.gap = fill.gap),
                strand(x),
                seqlengths = seqlengths(x))
    }
)

### FIXME: This is currently not doing a "punion" at all. It just appends
### the ranges in 'y' to their corresponding element in 'x'.
### 2 proposals for a more punion-like semantic:
###   (a) for (i in seq_len(length(x)))
###         x[[i]] <- punion(x[[i]], y[rep.int(i, length(x[[i]]))])
###   (b) for (i in seq_len(length(x)))
###         x[[i]] <- union(x[[i]], y[i])
### Note that behaviour (b) could also be considered a valid candidate for
### a union,GRangesList,GRanges method (which we don't have at the moment).
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
        ## TODO: 'seqinfo(x)' and 'seqinfo(y)' need to be merged. See TODO
        ## in "intersect" method at the top of this file for more details.
        ## Then the setequal check on the seqnames levels above won't be
        ## necessary anymore.
        GRanges(seqnames(x), ansRanges, ansStrand, seqlengths = seqlengths(x))
    }
)

### TODO: Like for "punion", the semantic of the
### pintersect,GRangesList,GRanges method should simply derive from
### the semantic of the pintersect,GRanges,GRanges. Then both, the
### implementation and documentation will be much easier to understand.
### It's hard to guess what's the current semantic of this method is by
### just looking at the code below, but it doesn't seem to be one of:
###   (a) for (i in seq_len(length(x)))
###         x[[i]] <- pintersect(x[[i]], y[rep.int(i, length(x[[i]]))])
###   (b) for (i in seq_len(length(x)))
###         x[[i]] <- intersect(x[[i]], y[i])
### It seems to be close to (b) but with special treatment of the "*"
### strand value in 'y'.
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
        ## TODO: 'seqinfo(x)' and 'seqinfo(y)' need to be merged. See TODO
        ## in "intersect" method at the top of this file for more details.
        ## Then the setequal check on the seqnames levels above won't be
        ## necessary anymore.
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
        ans <- GRanges(seqnames(x), ansRanges, ansStrand)
        ## TODO: 'seqinfo(x)' and 'seqinfo(y)' need to be merged. See TODO
        ## in "intersect" method at the top of this file for more details.
        ## Then the setequal check on the seqnames levels above won't be
        ## necessary anymore.
        ans@seqinfo <- x@seqinfo
        ans
    }
)

### TODO: Review the semantic of this method (see previous TODO's for
### "punion" and "pintersect" methods for GRanges,GRangesList).
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
          GRanges(ansSeqnames, unlist(ansRanges, use.names = FALSE), ansStrand)
        ## TODO: 'seqinfo(x)' and 'seqinfo(y)' need to be merged. See TODO
        ## in "intersect" method at the top of this file for more details.
        ## Then the setequal check on the seqnames levels above won't be
        ## necessary anymore.
        ansGRanges@seqinfo <- x@seqinfo
        new2("GRangesList",
             elementMetadata = new("DataFrame", nrows = length(x)),
             unlistData = ansGRanges, partitioning = ansRanges@partitioning,
             check=FALSE)
    }
)

