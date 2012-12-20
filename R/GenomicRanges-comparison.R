### =========================================================================
### Comparing and ordering genomic ranges
### -------------------------------------------------------------------------
###
### I. UNIQUE AND DUPLICATED ELEMENTS WITHIN A GenomicRanges OBJECT
### ---------------------------------------------------------------
### Two elements of a GenomicRanges object (i.e. two genomic ranges) are
### considered equal iff they are on the same underlying sequence and strand,
### and have the same start and width. duplicated() and unique() on a
### GenomicRanges object are conforming to this.
###
### II. ORDERING THE ELEMENTS OF A GenomicRanges OBJECT
### ---------------------------------------------------
### The "natural order" for the elements of a GenomicRanges object is to order
### them (a) first by sequence level, (b) then by strand, (c) then by start,
### (d) and finally by width.
### This way, the space of genomic ranges is totally ordered.
### Note that the "reduce" method for GenomicRanges uses this "natural order"
### implicitly. Also, note that, because we already do (c) and (d) for regular
### ranges, genomic ranges that belong to the same underlying sequence
### and strand are ordered like regular ranges.
### order(), sort(), and rank() on a GenomicRanges object are using this
### "natural order".
###
### III. ELEMENT-WISE (AKA "PARALLEL") COMPARISON OF 2 GenomicRanges OBJECTS
### ------------------------------------------------------------------------
### We want the "==", "!=", "<=", ">=", "<" and ">" operators between 2
### GenomicRanges objects to be compatible with the "natural order" defined
### previously. Defining those operators when the 2 objects have *identical*
### seqlevels() is straighforward but we can in fact extend this comparison
### to the following situation:
###   (A) 'e1' and 'e2' have compatible sets of underlying sequences, that is,
###       'seqinfo(e1)' and 'seqinfo(e2)' can be merged.
###   (B) 'seqlevels(e1)' and 'seqlevels(e2)' are in the same order. Note that
###       (A) guarantees that the seqlevels of one is a subset of the seqlevels
###       of the other. (B) is saying that this subset should be a subsequence.
### Pre-comparison step: if (A) and (B) are satisfied, then the 2 seqinfo() are
### merged and the seqlevels() of the result is assigned back to each object
### to compare. This is a way to have 2 objects with identical seqlevels()
### before the comparison can actually be performed and meaningful.
### The reason (B) is required for the pre-comparison step is because we want
### this step to preserve the original order of the seqlevels() in *both*
### objects. Without this precaution, the expected anti-symetric property of
### some operators would not be satisfied e.g. 'any(e1 < e2 & e2 < e1)' could
### be TRUE.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### compare()
###
### Doing 'compare(x, y)' on 2 vector-like objects 'x' and 'y' of length 1
### must return an integer less than, equal to, or greater than zero if the
### single element in 'x' is considered to be respectively less than, equal
### to, or greater than the single element in 'y'.
###
### On 2 GenomicRanges objects, it returns one of the 13 predefined codes
### (>= -6 and <= 6) used by the method for Ranges objects when 'x[i]' and
### 'y[i]' are on the same space (i.e. on the same underlying sequence and
### strand). Otherwise, it returns a code that is < -6 if 'x[i] < y[i]',
### and > 6 if 'x[i] > y[i]'.
### See '?compare' (in IRanges) for the 13 predefined codes.
.GenomicRanges.compare <- function(x, y)
{
    ## Pre-comparison step (see above for details).
    ## merge() will fail if 'x' and 'y' don't have compatible underlying
    ## sequences.
    seqinfo <- merge(seqinfo(x), seqinfo(y))
    seqlevels <- seqlevels(seqinfo)
    if (any(diff(match(seqlevels(y), seqlevels)) < 0L))
        stop("the 2 objects to compare have seqlevels in incompatible orders")
    ## This should only insert new seqlevels in the existing ones i.e. it
    ## should NEVER drop or reorder existing levels
    seqlevels(x) <- seqlevels(y) <- seqlevels
    a <- compare(ranges(x), ranges(y))
    b <- as.integer(strand(x)) - as.integer(strand(y))
    c <- as.integer(seqnames(x)) - as.integer(seqnames(y))
    ## Note that sign() always returns a numeric vector, even on an integer
    ## vector.
    a + 13L * as.integer(sign(b) + 3L * sign(c))
}

### The "compare" method.
setMethod("compare", c("GenomicRanges", "GenomicRanges"),
    function(x, y) .GenomicRanges.compare(x, y)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Element-wise (aka "parallel") comparison of 2 GenomicRanges objects.
###
### We only need to implement "==" and "<=" methods. The other comparison
### binary operators (!=, >=, <, >) will then work out-of-the-box on
### GenomicRanges objects thanks to the methods for Vector objects.
###

setMethod("==", signature(e1="GenomicRanges", e2="GenomicRanges"),
    function(e1, e2) { .GenomicRanges.compare(e1, e2) == 0L }
)

setMethod("<=", signature(e1="GenomicRanges", e2="GenomicRanges"),
    function(e1, e2) { .GenomicRanges.compare(e1, e2) <= 0L }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### duplicated()
###
### unique() will work out-of-the-box on a GenomicRanges object thanks to the
### method for Vector objects.
###

.duplicated.GenomicRanges <- function(x, incomparables=FALSE, fromLast=FALSE,
                                      method=c("auto", "quick", "hash"))
{
    if (!identical(incomparables, FALSE))
        stop("\"duplicated\" method for GenomicRanges objects ",
             "only accepts 'incomparables=FALSE'")
    IRanges:::duplicatedIntegerQuads(as.factor(seqnames(x)),
                                     as.factor(strand(x)),
                                     start(x), width(x),
                                     fromLast=fromLast, method=method)
}
### S3/S4 combo for duplicated.GenomicRanges
duplicated.GenomicRanges <- function(x, incomparables=FALSE, ...)
    .duplicated.GenomicRanges(x, incomparables=incomparables, ...)
setMethod("duplicated", "GenomicRanges", .duplicated.GenomicRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### match()
###
### %in% will work out-of-the-box on GenomicRanges objects thanks to the
### method for Vector objects.
###

### Unfortunately, the early version of this method was doing overlaps, not
### equality. We temporarily add the 'match.if.overlap' argument so the old
### behavior is still available.
### TODO: Deprecate 'match.if.overlap' arg in BioC 2.13.
### TODO: Remove 'match.if.overlap' arg in BioC 2.14.
setMethod("match", c("GenomicRanges", "GenomicRanges"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL,
                       method=c("auto", "quick", "hash"),
                       match.if.overlap=FALSE)
    {
        if (!isSingleNumberOrNA(nomatch))
            stop("'nomatch' must be a single number or NA")
        if (!is.integer(nomatch))
            nomatch <- as.integer(nomatch)
        if (!is.null(incomparables))
            stop("\"match\" method for GenomicRanges objects ",
                 "only accepts 'incomparables=NULL'")
        if (missing(match.if.overlap))
            warning(IRanges:::match.if.overlap.warning.msg("GenomicRanges"))
        if (!isTRUEorFALSE(match.if.overlap))
            stop("'match.if.overlap' must be TRUE or FALSE")
        if (match.if.overlap) {
            ans <- findOverlaps(x, table, select="first")
            if (!is.na(nomatch) && anyMissing(ans))
                ans[is.na(ans)] <- nomatch
            return(ans)
        }
        ## Equivalent to (but faster than):
        ##     findOverlaps(x, table, type="equal", select="first")
        ## except when 'x' and 'table' both contain empty ranges.
        IRanges:::matchIntegerQuads(as.factor(seqnames(x)),
                                    as.factor(strand(x)),
                                    start(x), width(x),
                                    as.factor(seqnames(table)),
                                    as.factor(strand(table)),
                                    start(table), width(table),
                                    nomatch=nomatch, method=method)
    }
)

### The only reason for overriding the method for Vector objects is to issue
### the warning.
### TODO: Remove this method in BioC 2.14 when the 'match.if.overlap' arg of
### match() is gone.
setMethod("%in%", c("GenomicRanges", "GenomicRanges"),
    function(x, table)
    {
        warning(IRanges:::`%in%.warning.msg`("GenomicRanges"))
        !is.na(match(x, table, match.if.overlap=FALSE))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### order() and related methods.
###
### The "order" and "rank" methods for GenomicRanges objects are consistent
### with the order implied by compare().
### sort() will work out-of-the-box on a GenomicRanges object thanks to the
### method for Vector objects.
###

setMethod("order", "GenomicRanges",
    function(..., na.last=TRUE, decreasing=FALSE)
    {
        if (!isTRUEorFALSE(decreasing))
            stop("'decreasing' must be TRUE or FALSE")
        ## All arguments in '...' are guaranteed to be GenomicRanges objects.
        args <- list(...)
        if (length(args) == 1L) {
            x <- args[[1L]]
            return(IRanges:::orderIntegerQuads(as.factor(seqnames(x)),
                                               as.factor(strand(x)),
                                               start(x), width(x),
                                               decreasing=decreasing))
        }
        order_args <- vector("list", 4L*length(args))
        idx <- 4L*seq_len(length(args))
        order_args[idx - 3L] <- lapply(args,
                                       function(x) as.factor(seqnames(x)))
        order_args[idx - 2L] <- lapply(args,
                                       function(x) as.factor(strand(x)))
        order_args[idx - 1L] <- lapply(args, start)
        order_args[idx] <- lapply(args, width)
        do.call(order, c(order_args,
                         list(na.last=na.last, decreasing=decreasing)))
    }
)

setMethod("rank", "GenomicRanges",
    function(x, na.last=TRUE,
             ties.method=c("average", "first", "random", "max", "min"))
    {
        if (!missing(ties.method) && !identical(ties.method, "first"))
            stop("only 'ties.method=\"first\"' is supported ",
                 "when ranking genomic ranges")
        oo <- order(x)
        ## 'ans' is the reverse permutation of 'oo'
        ans <- integer(length(oo))
        ans[oo] <- seq_len(length(oo))
        ans
    }
)

