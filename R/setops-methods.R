### =========================================================================
### Set operations
### -------------------------------------------------------------------------

### TODO: What's the impact of circularity on the set operations?


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 3 low-level helper functions.
###

### A fast implementation of 'mendoapply(FUN, x, y)' for GRangesList objects.
### Assume 'x' and 'y' are 2 GRangesList objects (not checked) of the same
### length (checked).
.fast_binary_mendoapply <- function(FUN, x, y, ...)
{
    FUN <- match.fun(FUN)
    if (length(x) != length(y))
        stop("'x' and 'y' must have the same length")
    seqinfo(x) <- merge(seqinfo(x), seqinfo(y))
    seqlevels(y) <- seqlevels(x)
    xgr <- deconstructGRLintoGR(x)
    ygr <- deconstructGRLintoGR(y)
    seqinfo(xgr) <- suppressWarnings(merge(seqinfo(xgr), seqinfo(ygr)))
    seqlevels(ygr) <- seqlevels(xgr)
    gr <- FUN(xgr, ygr, ...)
    reconstructGRLfromGR(gr, x)
}

### Both return a named integer vector where the names are guaranteed to be
### 'seqlevels(x)'.
###

.minStartPerGRangesSequence <- function(x)
{
    cil <- splitAsList(start(x), seqnames(x))  # CompressedIntegerList object
    ## The 4 lines below are equivalent to:
    ##   ans <- min(cil)
    ##   ans[elementNROWS(v) == 0L] <- NA_integer_
    ## but much faster!
    ## TODO: Replace with the above, but only when the "min" method for
    ## CompressedIntegerList objects (implemented in the IRanges package)
    ## is as fast as the "viewMins" method for XIntegerViews objects
    ## (implemented in C in the XVector package). Ideally, the 2 methods
    ## should share the same underlying code.
    v <- Views(cil@unlistData, cil@partitioning)  # XIntegerViews object
    ans <- viewMins(v)
    ans[width(v) == 0L] <- NA_integer_
    names(ans) <- names(v)
    ans
}

.maxEndPerGRangesSequence <- function(x)
{
    cil <- splitAsList(end(x), seqnames(x))  # CompressedIntegerList object
    ## The 4 lines below are equivalent to:
    ##   ans <- max(cil)
    ##   ans[elementNROWS(v) == 0L] <- NA_integer_
    ## but much faster!
    ## TODO: Replace with the above, but only when the "max" method for
    ## CompressedIntegerList objects (implemented in the IRanges package)
    ## is as fast as the "viewMaxs" method for XIntegerViews objects
    ## (implemented in C in the XVector package). Ideally, the 2 methods
    ## should share the same underlying code.
    v <- Views(cil@unlistData, cil@partitioning)  # XIntegerViews object
    ans <- viewMaxs(v)
    ans[width(v) == 0L] <- NA_integer_
    names(ans) <- names(v)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### union(), intersect(), setdiff()
###

### Always return a GRanges *instance* whatever GenomicRanges derivatives are
### passed to it (e.g. GPos or GNCList), so does NOT act like an endomorphism
### in general.
setMethod("union", c("GenomicRanges", "GenomicRanges"),
    function(x, y, ignore.strand=FALSE)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        x <- granges(x)
        y <- granges(y)
        if (ignore.strand)
            strand(x) <- strand(y) <- "*"
        reduce(c(x, y), drop.empty.ranges=TRUE)
    }
)

.unsupported_union <- function(x, y)
    stop("union() between a ", class(x), " and a ", class(y), " object ",
         "is not supported")

setMethod("union", c("GenomicRanges", "Vector"), .unsupported_union)
setMethod("union", c("Vector", "GenomicRanges"), .unsupported_union)

setMethod("union", c("GRangesList", "GRangesList"),
    function(x, y, ...) .fast_binary_mendoapply("union", x, y, ...)
)

### Always return a GRanges *instance* whatever GenomicRanges derivatives are
### passed to it (e.g. GPos or GNCList), so does NOT act like an endomorphism
### in general.
setMethod("intersect", c("GenomicRanges", "GenomicRanges"),
    function(x, y, ignore.strand=FALSE)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        x <- granges(x)
        y <- granges(y)
        if (ignore.strand)
            strand(x) <- strand(y) <- "*"
        seqinfo(x) <- merge(seqinfo(x), seqinfo(y))
        ## If merge() is going to issue a warning, we don't want to get
        ## it twice.
        seqinfo(y) <- suppressWarnings(merge(seqinfo(y), seqinfo(x)))
        seqlengths <- seqlengths(x)
        ## If the length of a sequence is unknown (NA), then we use
        ## the max end value found on that sequence in 'x' or 'y'.
        seqlengths[is.na(seqlengths)] <-
            .maxEndPerGRangesSequence(c(x, y))[is.na(seqlengths)]
        setdiff(x, gaps(y, end=seqlengths))
    }
)

.unsupported_intersect <- function(x, y)
    stop("intersect() between a ", class(x), " and a ", class(y), " object ",
         "is not supported")

setMethod("intersect", c("GenomicRanges", "Vector"), .unsupported_intersect)
setMethod("intersect", c("Vector", "GenomicRanges"), .unsupported_intersect)

setMethod("intersect", c("GRangesList", "GRangesList"),
    function(x, y, ...) .fast_binary_mendoapply("intersect", x, y, ...)
)

### Always return a GRanges *instance* whatever GenomicRanges derivatives are
### passed to it (e.g. GPos or GNCList), so does NOT act like an endomorphism
### in general.
setMethod("setdiff", c("GenomicRanges", "GenomicRanges"),
    function(x, y, ignore.strand=FALSE)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        x <- granges(x)
        y <- granges(y)
        if (ignore.strand)
            strand(x) <- strand(y) <- "*"
        seqinfo(x) <- merge(seqinfo(x), seqinfo(y))
        ## If merge() is going to issue a warning, we don't want to get
        ## it twice.
        seqinfo(y) <- suppressWarnings(merge(seqinfo(y), seqinfo(x)))
        seqlengths <- seqlengths(x)
        ## If the length of a sequence is unknown (NA), then we use
        ## the max end value found on that sequence in 'x' or 'y'.
        seqlengths[is.na(seqlengths)] <-
            .maxEndPerGRangesSequence(c(x, y))[is.na(seqlengths)]
        gaps(union(gaps(x, end=seqlengths), y), end=seqlengths)
    }
)

.unsupported_setdiff <- function(x, y)
    stop("setdiff() between a ", class(x), " and a ", class(y), " object ",
         "is not supported")

setMethod("setdiff", c("GenomicRanges", "Vector"), .unsupported_setdiff)
setMethod("setdiff", c("Vector", "GenomicRanges"), .unsupported_setdiff)

setMethod("setdiff", c("GRangesList", "GRangesList"),
    function(x, y, ...) .fast_binary_mendoapply("setdiff", x, y, ...)
)


### =========================================================================
### Parallel set operations
### -------------------------------------------------------------------------

## check for equality without requiring identical level sets
## instead, one level set must be subset of the other, like merge,Seqinfo
compatibleSeqnames <- function(x, y) {
  if (length(x) != length(y))
    stop("'x' and 'y' must be of equal length")
  if (!is(x, "Rle") || !is(y, "Rle"))
    stop("'x' and 'y' must be Rle objects")
  xLevels <- levels(x)
  yLevels <- levels(y)
  if (length(xLevels) > length(yLevels))
    diffLevels <- setdiff(yLevels, xLevels)
  else diffLevels <- setdiff(xLevels, yLevels)
  if (length(diffLevels))
    stop("Level set of 'x' must be subset of that of 'y', or vice versa")
  runValue(x) <- as.character(runValue(x))
  runValue(y) <- as.character(runValue(y))
  x == y
}

allCompatibleSeqnamesAndStrand <- function(x, y) {
  all(compatibleSeqnames(seqnames(x), seqnames(y)) &
      compatibleStrand(strand(x), strand(y)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### punion()
###

setMethod("punion", c("GRanges", "GRanges"),
    function(x, y, fill.gap=FALSE, ignore.strand=FALSE)
    {
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand) 
            strand(y) <- strand(x)
        mcols(x) <- NULL
        seqinfo(x) <- merge(seqinfo(x), seqinfo(y))
        if (!allCompatibleSeqnamesAndStrand(x, y))
            stop("'x' and 'y' elements must have compatible 'seqnames' ",
                 "and 'strand' values")
        ranges(x) <- punion(ranges(x), ranges(y), fill.gap=fill.gap)
        x
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
    function(x, y, fill.gap=FALSE)
    {
        n <- length(x)
        if (n != length(y)) 
            stop("'x' and 'y' must have the same length")
        mcols(x@unlistData) <- NULL
        mcols(y) <- NULL
        ans <-
          split(c(x@unlistData, y), 
                c(Rle(seq_len(n), elementNROWS(x)), Rle(seq_len(n))))
        names(ans) <- names(x)
        ans
    }
)

setMethod("punion", c("GRanges", "GRangesList"),
    function(x, y, ...) callGeneric(y, x, ...)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pintersect()
###

### 'x' and 'y' must have the same length. An 'y' of length 1 is accepted and
### is recycled to the length of 'x'.
### Return a GRanges object. If 'drop.nohit.ranges' is FALSE (the default) the
### returned object is parallel to 'x' with the original ranges modified
### (seqnames, names, and metadata are untouched), so this can be seen as an
### "intra-range transformation". The original metadata columns are propagated
### and a new "hit" column added to them to indicate elements in 'x' that
### intersect with the corresponding element in 'y'. For these elements the
### range in the returned object is guaranteed to be a subrange of the
### original ranges. If 'ignore.strand' or 'strict.strand' is TRUE, then the
### returned GRanges has the same strand as 'x'. Otherwise, for elements in 'x'
### that are on the "*" strand and hit the corresponding element in 'y', it
### has the strand of 'y'.
### If 'drop.nohit.ranges' is TRUE everything is the same except that elements
### in 'x' that don't intersect with the corresponding element in 'y' are
### removed from the result (so the result is no more parallel to the input)
### and no "hit" metadata column is added to it.
.pintersect_GRanges_GRanges <- function(x, y, drop.nohit.ranges=FALSE,
                                              ignore.strand=FALSE,
                                              strict.strand=FALSE)
{
    stopifnot(is(x, "GRanges"), is(y, "GRanges"))
    if (length(y) != length(x)) {
        if (length(y) != 1L)
            stop(wmsg("'y' must have the length of 'x' or length 1"))
        y <- rep.int(y, length(x))
    }
    if (!isTRUEorFALSE(drop.nohit.ranges))
        stop(wmsg("'drop.nohit.ranges' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(ignore.strand))
        stop(wmsg("'ignore.strand' must be TRUE or FALSE"))
    if (!isTRUEorFALSE(strict.strand))
        stop(wmsg("'strict.strand' must be TRUE or FALSE"))
    if (ignore.strand) {
        if (strict.strand)
            warning(wmsg("'strict.strand' is ignored ",
                         "when 'ignore.strand' is TRUE"))
        strand(y) <- strand(x)
    } else if (!strict.strand) {
        x_strand <- strand(x)
        idx <- strand(x) == "*" & strand(y) != "*"
        strand(x)[idx] <- strand(y)[idx]
        idx <- strand(y) == "*" & strand(x) != "*"
        strand(y)[idx] <- strand(x)[idx]
    }

    x_seqlevels <- seqlevels(x)
    ## Check compatibility of underlying genomes.
    seqinfo(x) <- merge(seqinfo(x), seqinfo(y))
    ## Align the seqlevels of 'y' to those of 'x' so that comparing the
    ## seqnames of the 2 objects thru as.integer(seqnames()) is meaningful).
    seqlevels(y) <- seqlevels(x)
    ## Restore original seqlevels on 'x' (this preserves their order).
    seqlevels(x) <- x_seqlevels

    new_start <- pmax.int(start(x), start(y))
    new_end <- pmin.int(end(x), end(y))

    is_hit <- as.integer(seqnames(x)) == as.integer(seqnames(y)) &
              strand(x) == strand(y) & new_end >= new_start - 1L

    if (drop.nohit.ranges) {
        x <- extractROWS(x, is_hit)
        new_start <- extractROWS(new_start, is_hit)
        new_end <- extractROWS(new_end, is_hit)
        new_names <- names(x)
    } else {
        ## For elements in 'x' and 'y' that don't hit each other, we return
        ## an artificial zero-width intersection that starts on 'start(x)'.
        nohit_idx <- which(!is_hit)
        nohit_start <- start(x)[nohit_idx]
        new_start[nohit_idx] <- nohit_start
        new_end[nohit_idx] <- nohit_start - 1L
        new_names <- names(x)
        mcols(x)$hit <- as.logical(is_hit)
        if (!(ignore.strand || strict.strand))
            strand(x)[nohit_idx] <- x_strand[nohit_idx]
    }
    ranges(x) <- IRanges(new_start, new_end, names=new_names)
    x
}

setMethod("pintersect", c("GRanges", "GRanges"),
    function(x, y, drop.nohit.ranges=FALSE,
                   ignore.strand=FALSE, strict.strand=FALSE)
        .pintersect_GRanges_GRanges(x, y, drop.nohit.ranges=drop.nohit.ranges,
                                          ignore.strand=ignore.strand,
                                          strict.strand=strict.strand)
)

### This is equivalent to 'mendoapply(pintersect, x, y)'. Therefore the
### returned GRangesList object has the same shape as 'x'.
### To get the 'mendoapply(intersect, x, y)' behavior, the user should use
### 'strict.strand=TRUE' and call 'reduce( , drop.empty.ranges=TRUE)' on
### the returned object.
.pintersect_GRangesList_GRanges <- function(x, y, drop.nohit.ranges=FALSE,
                                                  ignore.strand=FALSE,
                                                  strict.strand=FALSE)
{
    stopifnot(is(x, "GRangesList"), is(y, "GRanges"))
    if (length(y) != length(x)) {
        if (length(y) != 1L)
            stop(wmsg("'y' must have the length of 'x' or length 1"))
        y <- rep.int(y, length(x))
    }
    if (!isTRUEorFALSE(drop.nohit.ranges))
        stop(wmsg("'drop.nohit.ranges' must be TRUE or FALSE"))
    unlisted_x <- unlist(x, use.names=FALSE)
    y2 <- rep.int(y, elementNROWS(x))
    ## 'unlisted_ans' parallel to 'x'.
    unlisted_ans <- .pintersect_GRanges_GRanges(unlisted_x, y2,
                                        drop.nohit.ranges=FALSE,
                                        ignore.strand=ignore.strand,
                                        strict.strand=strict.strand)
    if (drop.nohit.ranges) {
        is_hit <- relist(mcols(unlisted_ans, use.names=FALSE)$hit, x)
        mcols(unlisted_ans)$hit <- NULL
        ans <- relist(unlisted_ans, x)[is_hit]
    } else {
        ans <- relist(unlisted_ans, x)
    }
    ans
}

setMethod("pintersect", c("GRangesList", "GRanges"),
    function(x, y, drop.nohit.ranges=FALSE,
                   ignore.strand=FALSE, strict.strand=FALSE)
        .pintersect_GRangesList_GRanges(x, y,
                                        drop.nohit.ranges=drop.nohit.ranges,
                                        ignore.strand=ignore.strand,
                                        strict.strand=strict.strand)
)

setMethod("pintersect", c("GRanges", "GRangesList"),
    function(x, y, drop.nohit.ranges=FALSE,
                   ignore.strand=FALSE, strict.strand=FALSE)
        .pintersect_GRangesList_GRanges(y, x,
                                        drop.nohit.ranges=drop.nohit.ranges,
                                        ignore.strand=ignore.strand,
                                        strict.strand=strict.strand)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### psetdiff()
###

setMethod("psetdiff", c("GRanges", "GRanges"),
    function(x, y, ignore.strand=FALSE)
    {
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand)
            strand(y) <- strand(x)
        mcols(x) <- NULL
        seqinfo(x) <- merge(seqinfo(x), seqinfo(y))
        ok <- compatibleSeqnames(seqnames(x), seqnames(y)) &
              compatibleStrand(strand(x), strand(y))
        ## Update the ranges.
        ansRanges <- ranges(x)
        ansRanges <- replaceROWS(ansRanges, ok,
                         callGeneric(extractROWS(ranges(x), ok),
                                     extractROWS(ranges(y), ok)))
        ranges(x) <- ansRanges
        ## Update the strand.
        ansStrand <- strand(x)
        resolveStrand <- as(ansStrand == "*", "IRanges")
        if (length(resolveStrand) > 0L)
            ansStrand[IRanges:::unlist_as_integer(resolveStrand)] <-
              extractROWS(strand(y), resolveStrand)
        strand(x) <- ansStrand
        x
    }
)

### TODO: Review the semantic of this method (see previous TODO's for
### "punion" and "pintersect" methods for GRanges,GRangesList).
setMethod("psetdiff", c("GRanges", "GRangesList"),
    function(x, y, ignore.strand=FALSE)
    {
        ansSeqinfo <- merge(seqinfo(x), seqinfo(y))
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        ok <- compatibleSeqnames(rep(seqnames(x), elementNROWS(y)),
                                 seqnames(y@unlistData))
        if (!ignore.strand)
          ok <- ok & compatibleStrand(rep(strand(x), elementNROWS(y)),
                                      strand(y@unlistData))
        if (!all(ok)) {
            ok <-
              new2("CompressedLogicalList", unlistData=as.vector(ok),
                   partitioning=y@partitioning)
            y <- y[ok]
        }
        ansRanges <- gaps(ranges(y), start=start(x), end=end(x))
        ansSeqnames <- rep(seqnames(x), elementNROWS(ansRanges))
        ansStrand <- rep(strand(x), elementNROWS(ansRanges))
        ansGRanges <-
          GRanges(ansSeqnames, unlist(ansRanges, use.names=FALSE), ansStrand)
        seqinfo(ansGRanges) <- ansSeqinfo
        relist(ansGRanges, PartitioningByEnd(ansRanges))
    }
)

setMethod("pgap", c("GRanges", "GRanges"),
    function(x, y, ignore.strand=FALSE, ...)
    {
        if (length(x) != length(y)) 
            stop("'x' and 'y' must have the same length")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand) 
            strand(y) <- strand(x) 
        mcols(x) <- NULL
        seqinfo(x) <- merge(seqinfo(x), seqinfo(y))
        if (!allCompatibleSeqnamesAndStrand(x, y))
            stop("'x' and 'y' elements must have compatible 'seqnames' ",
                 "and 'strand' values")
        ranges(x) <- callGeneric(ranges(x), ranges(y), ...)
        x
    }
)

