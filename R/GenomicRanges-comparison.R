### =========================================================================
### Ordering and comparing genomic ranges
### -------------------------------------------------------------------------
###
### I. UNIQUE AND DUPLICATED ELEMENTS WITHIN A GenomicRanges OBJECT
### ---------------------------------------------------------------
### Two elements of a GenomicRanges object (i.e. two genomic ranges) are
### considered equal iff they are on the same underlying sequence and strand,
### and have the same start and width.
### The "duplicated" and "unique" methods for GenomicRanges objects are
### using this equality.
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
### The "order", "sort" and "rank" methods for GenomicRanges objects are using
### this "natural order".
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
### I. UNIQUE AND DUPLICATED ELEMENTS WITHIN A GenomicRanges OBJECT.
###

setMethod("duplicated", "GenomicRanges",
          function(x, incomparables=FALSE, fromLast=FALSE,
                   method=c("auto", "quick", "hash"), ...)
          {
            if (!identical(incomparables, FALSE))
              stop("\"duplicated\" method for GenomicRanges objects ",
                   "only accepts 'incomparables=FALSE'")
            IRanges:::duplicatedIntegerQuads(as.factor(seqnames(x)),
                                             as.factor(strand(x)),
                                             start(x), width(x),
                                             fromLast=fromLast, method=method)
          }
          )


### Relies on a "[" method for 'x'.
setMethod("unique", "GenomicRanges",
    function(x, incomparables=FALSE, fromLast=FALSE,
             method=c("auto", "quick", "hash"), ...)
    {
        x[!duplicated(x, incomparables=incomparables,
                         fromLast=fromLast, method=method, ...)]
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### II. ORDERING THE ELEMENTS OF A GenomicRanges OBJECT
###

setMethod("order", "GenomicRanges",
    function(..., na.last=TRUE, decreasing=FALSE)
    {
        if (!isTRUEorFALSE(decreasing))
            stop("'decreasing' must be TRUE or FALSE")
        ## all arguments in '...' are guaranteed to be GenomicRanges objects
        args <- list(...)
        if (length(args) == 1L) {
            x <- args[[1L]]
            return(IRanges:::orderIntegerQuads(as.factor(seqnames(x)),
                                               as.factor(strand(x)),
                                               start(x), width(x),
                                               decreasing = decreasing))
        }
        order_args <- vector("list", 4L*length(args))
        idx <- 4L*seq_len(length(args))
        order_args[idx - 3L] <- lapply(args,
                                       function(x) as.factor(seqnames(x)))
        order_args[idx - 2L] <- lapply(args,
                                       function(x) as.factor(strand(x)))
        order_args[idx - 1L] <- lapply(args, start)
        order_args[idx] <- lapply(args, width)
        do.call(order, c(order_args, list(decreasing=decreasing)))
    }
)

### Relies on a "[" method for 'x'.
### TODO: Define a "sort" method for ANY (in IRanges) instead of this.
setMethod("sort", "GenomicRanges",
    function(x, decreasing=FALSE, ...) 
    {
        x[order(x, decreasing=decreasing)]
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### III. ELEMENT-WISE (AKA "PARALLEL") COMPARISON OF 2 GenomicRanges OBJECTS.
###

### This is the core function for parallel comparison of 2 GenomicRanges
### objects. It returns a numeric vector of the length of the longest
### GenomicRanges where the i-th value is less than, equal to, or greater than
### zero if the i-th element in 'e1' is respectively less than, equal to, or
### greater than the i-th element in 'e2'.
.GenomicRanges.pcompar <- function(e1, e2)
{
    ## Pre-comparison step (see above for details).
    ## merge() will fail if 'e1' and 'e2' don't have compatible
    ## underlying sequences.
    seqinfo <- merge(seqinfo(e1), seqinfo(e2))
    seqlevels <- seqlevels(seqinfo)
    if (any(diff(match(seqlevels(e2), seqlevels)) < 0L))
        stop("the 2 objects to compare have seqlevels in incompatible orders")
    ## This should only insert new seqlevels in the existing ones i.e. it
    ## should NEVER drop or reorder existing levels
    seqlevels(e1) <- seqlevels(e2) <- seqlevels(seqinfo)
    a <- as.integer(seqnames(e1)) - as.integer(seqnames(e2))
    b <- as.integer(strand(e1)) - as.integer(strand(e2))
    c <- start(e1) - start(e2)
    d <- width(e1) - width(e2)
    ## Note that sign() always returns a numeric vector, even on an integer
    ## vector.
    as.integer(8L * sign(a) + 4L * sign(b) + 2L * sign(c) + sign(d))
}

### There is a "!=" method for ANY,ANY defined in IRanges.
setMethod("==", signature(e1="GenomicRanges", e2="GenomicRanges"),
    function(e1, e2) { .GenomicRanges.pcompar(e1, e2) == 0L }
)

setMethod("<=", signature(e1="GenomicRanges", e2="GenomicRanges"),
    function(e1, e2) { .GenomicRanges.pcompar(e1, e2) <= 0L }
)

### TODO: Define methods for ANY,ANY (in IRanges) instead of the 3 following
### methods.
setMethod(">=", signature(e1="GenomicRanges", e2="GenomicRanges"),
    function(e1, e2) { e2 <= e1 }
)
setMethod("<", signature(e1="GenomicRanges", e2="GenomicRanges"),
    function(e1, e2) { !(e2 <= e1) }
)
setMethod(">", signature(e1="GenomicRanges", e2="GenomicRanges"),
    function(e1, e2) { !(e1 <= e2) }
)

