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
### is.unsorted(), order(), sort(), and rank() on a GenomicRanges object
### adhere to this "natural order".
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


setMethod("pcompareRecursively", "GenomicRanges", function(x) FALSE)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pcompare()
###
### Doing 'pcompare(x, y)' on 2 vector-like objects 'x' and 'y' of length 1
### must return an integer less than, equal to, or greater than zero if the
### single element in 'x' is considered to be respectively less than, equal
### to, or greater than the single element in 'y'.
###
### On 2 GenomicRanges objects, it returns one of the 13 predefined codes
### (>= -6 and <= 6) used by the method for IntegerRanges objects when 'x[i]'
### and 'y[i]' are on the same space (i.e. on the same underlying sequence
### and strand). Otherwise, it returns a code that is < -6 if 'x[i] < y[i]',
### and > 6 if 'x[i] > y[i]'.
### See '?pcompare' (in IRanges) for the 13 predefined codes.
.pcompare_GenomicRanges <- function(x, y)
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
    a <- pcompare(ranges(x), ranges(y))
    b <- as.integer(strand(x)) - as.integer(strand(y))
    c <- as.integer(seqnames(x)) - as.integer(seqnames(y))
    ## Note that sign() always returns a numeric vector, even on an integer
    ## vector.
    a + 13L * as.integer(sign(b) + 3L * sign(c))
}

### The "pcompare" method.
setMethod("pcompare", c("GenomicRanges", "GenomicRanges"),
    function(x, y) .pcompare_GenomicRanges(x, y)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### duplicated()
###
### unique() will work out-of-the-box on a GenomicRanges object thanks to the
### method for Vector objects.
###

.duplicated.GenomicRanges <- function(x, incomparables=FALSE, fromLast=FALSE,
                                      nmax = NA,
                                      method=c("auto", "quick", "hash"))
{
    if (!identical(incomparables, FALSE))
        stop("\"duplicated\" method for GenomicRanges objects ",
             "only accepts 'incomparables=FALSE'")
    duplicatedIntegerQuads(as.factor(seqnames(x)), as.factor(strand(x)),
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

setMethod("match", c("GenomicRanges", "GenomicRanges"),
    function(x, table, nomatch=NA_integer_, incomparables=NULL,
                       method=c("auto", "quick", "hash"),
                       ignore.strand=FALSE)
    {
        if (!isSingleNumberOrNA(nomatch))
            stop("'nomatch' must be a single number or NA")
        if (!is.integer(nomatch))
            nomatch <- as.integer(nomatch)
        if (!is.null(incomparables))
            stop("\"match\" method for GenomicRanges objects ",
                 "only accepts 'incomparables=NULL'")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        ## Calling merge() is the way to check that 'x' and 'table' are based
        ## on the same reference genome.
        merge(seqinfo(x), seqinfo(table))
        x_seqnames <- relevelSeqnamesForMatch(x, table)
        if (ignore.strand) {
            x_strand <- integer(length(x))
            table_strand <- integer(length(table))
        } else {
            x_strand <- as.factor(strand(x))
            table_strand <- as.factor(strand(table))
        }
        ## Equivalent to (but faster than):
        ##     findOverlaps(x, table, type="equal", select="first")
        ## except when 'x' and 'table' both contain empty ranges.
        matchIntegerQuads(x_seqnames, x_strand,
                              start(x), width(x),
                          as.factor(seqnames(table)), table_strand,
                              start(table), width(table),
                          nomatch=nomatch, method=method)
    }
)

relevelSeqnamesForMatch <- function(x, table) {
  x_seqnames <- as.factor(seqnames(x))
  if (!hasHead(seqlevels(x), seqlevels(table)) &&
      !hasHead(seqlevels(table), seqlevels(x))) {
    x_seqnames <- factor(x_seqnames,
                         union(seqlevels(table), seqlevels(x)))
  }
  x_seqnames
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### selfmatch()
###

setMethod("selfmatch", "GenomicRanges",
    function(x, method=c("auto", "quick", "hash"), ignore.strand=FALSE)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        x_seqnames <- as.factor(seqnames(x))
        if (ignore.strand) {
            x_strand <- integer(length(x))
        } else {
            x_strand <- as.factor(strand(x))
        }
        selfmatchIntegerQuads(x_seqnames, x_strand, start(x), width(x),
                              method=method)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### order() and related methods.
###
### is.unsorted(), order(), sort(), rank() on GenomicRanges derivatives are
### consistent with the order implied by pcompare().
### is.unsorted() is a quick/cheap way of checking whether a GenomicRanges
### derivative is already sorted, e.g., called prior to a costly sort.
###

.GenomicRanges_as_integer_quads <- function(x, ignore.strand=FALSE)
{
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE of FALSE")
    a <- S4Vectors:::decodeRle(seqnames(x))
    if (ignore.strand) {
        b <- integer(length(x))
    } else {
        b <- S4Vectors:::decodeRle(strand(x))
    }
    c <- start(x)
    d <- width(x)
    list(a, b, c, d)
}

setMethod("is.unsorted", "GenomicRanges",
    function(x, na.rm=FALSE, strictly=FALSE, ignore.strand=FALSE)
    {
        if (!identical(na.rm, FALSE))
            warning("\"is.unsorted\" method for GenomicRanges objects ",
                    "ignores the 'na.rm' argument")
        if (!isTRUEorFALSE(strictly))
            stop("'strictly' must be TRUE of FALSE")
        ## It seems that creating the integer quads below is faster when
        ## 'x' is already sorted (TODO: Investigate why). Therefore, and
        ## somewhat counterintuitively, is.unsorted() can be faster when 'x'
        ## is already sorted (which, in theory, is the worst-case scenario
        ## because S4Vectors:::sortedIntegerQuads() will then need to take a
        ## full walk on 'x') than when it is unsorted (in which case
        ## S4Vectors:::sortedIntegerQuads() might stop walking on 'x' after
        ## checking its first 2 elements only -- the best-case scenario).
        quads <- .GenomicRanges_as_integer_quads(x, ignore.strand)
        !S4Vectors:::sortedIntegerQuads(quads[[1L]], quads[[2L]],
                                        quads[[3L]], quads[[4L]],
                                        strictly=strictly)
    }
)

### NOT exported but used in GenomicAlignments package.
order_GenomicRanges <- function(x, decreasing=FALSE, ignore.strand=FALSE)
{
    if (!isTRUEorFALSE(decreasing))
        stop("'decreasing' must be TRUE or FALSE")
    quads <- .GenomicRanges_as_integer_quads(x, ignore.strand)
    orderIntegerQuads(quads[[1L]], quads[[2L]],
                      quads[[3L]], quads[[4L]],
                      decreasing=decreasing)
}

### TODO: Support the 'ignore.strand' argument (the signature of the order()
### generic doesn't make this as straightforward as it could be).
### 'na.last' is pointless (GenomicRanges objects don't contain NAs) so is
### ignored.
### 'method' is also ignored at the moment.
setMethod("order", "GenomicRanges",
    function(..., na.last=TRUE, decreasing=FALSE,
                  method=c("auto", "shell", "radix"))
    {
        ## Turn off this warning for now since it triggers spurious warnings
        ## when calling sort() on a GRangesList object. The root of the
        ## problem is inconsistent defaults for 'na.last' between order() and
        ## sort(), as reported here:
        ##   https://stat.ethz.ch/pipermail/r-devel/2015-November/072012.html
        #if (!identical(na.last, TRUE))
        #    warning("\"order\" method for GenomicRanges objects ",
        #            "ignores the 'na.last' argument")
        if (!isTRUEorFALSE(decreasing))
            stop("'decreasing' must be TRUE or FALSE")
        ## All arguments in '...' are guaranteed to be GenomicRanges objects.
        args <- list(...)
        if (length(args) == 1L)
            return(order_GenomicRanges(args[[1L]], decreasing))
        order_args <- c(unlist(lapply(args, .GenomicRanges_as_integer_quads),
                               recursive=FALSE, use.names=FALSE),
                        list(na.last=na.last, decreasing=decreasing))
        do.call(order, order_args)
    }
)

### S3/S4 combo for sort.GenomicRanges
.sort.GenomicRanges <- function(x, decreasing=FALSE, ignore.strand=FALSE, by)
{
    if (missing(by)) {
        oo <- order_GenomicRanges(x, decreasing, ignore.strand)
    } else {
        if (!identical(ignore.strand, FALSE))
            warning("'ignore.strand' ignored when 'by' is specified")
        oo <- S4Vectors:::orderBy(by, x, decreasing=decreasing)
    }
    extractROWS(x, oo)
}
sort.GenomicRanges <- function(x, decreasing=FALSE, ...)
    .sort.GenomicRanges(x, decreasing=decreasing, ...)
setMethod("sort", "GenomicRanges", .sort.GenomicRanges)

setMethod("rank", "GenomicRanges",
    function(x, na.last=TRUE,
             ties.method=c("average", "first", "last", "random", "max", "min"),
             ignore.strand=FALSE)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE of FALSE")
        if (ignore.strand)
            x <- unstrand(x)
        callNextMethod(x, na.last=na.last, ties.method=ties.method)
    }
)

