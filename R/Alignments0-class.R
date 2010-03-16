### =========================================================================
### Alignments0 objects
### -------------------------------------------------------------------------
###

### First GappedAlignments implementation: Alignments0
### The genomic ranges of the alignments are stored in 3 slots: 'rname',
### 'strand' and 'rglist'.
### See http://wiki.fhcrc.org/bioc/Multiple_alignment_rep_v1 for the details
### of this class proposal.
setClass("Alignments0",
    contains="GappedAlignments",
    representation(
        rname="factor",               # character factor
        strand="raw",
        rglist="CompressedNormalIRangesList"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Some low-level helper functions for manipulating the @strand slot.
###

.logicalAsCompactRawVector <- function(x)
{
    if (!is.logical(x))
        stop("'x' must be a logical vector")
    .Call("logical_as_compact_raw_vector", x, PACKAGE="GenomicRanges")
}

.compactRawVectorAsLogical <- function(x, length.out)
{
    if (!is.raw(x))
        stop("'x' must be a raw vector")
    if (!isSingleNumber(length.out))
        stop("'length.out' must be a single number")
    if (!is.integer(length.out))
        length.out <- as.integer(length.out)
    .Call("compact_raw_vector_as_logical", x, length.out, PACKAGE="GenomicRanges")
}

.subsetCompactRawVector <- function(x, i)
{
    if (!is.raw(x))
        stop("'x' must be a raw vector")
    if (!is.integer(i))
        stop("'i' must be an integer vector")
    .Call("subset_compact_raw_vector", x, i, PACKAGE="GenomicRanges")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("rname", "Alignments0", function(x) x@rname)

setReplaceMethod("rname", "Alignments0",
    function(x, value)
    {
        x@rname <- normargRNameReplaceValue(x, value)
        x
    }
)

setMethod("strand", "Alignments0",
    function(x) strand(.compactRawVectorAsLogical(x@strand, length(x)))
)

setMethod("rglist", "Alignments0", function(x) x@rglist)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

Alignments0 <- function(rname=factor(), strand=NULL,
                        pos=integer(0), cigar=character(0))
{
    if (!is.factor(rname) || !is.character(levels(rname))) {
        if (!is.character(rname))
            stop("'rname' must be a character vector/factor")
        rname <- as.factor(rname)
    }
    if (any(is.na(rname)))
        stop("'rname' cannot have NAs")
    if (is.null(strand)) {
        if (length(rname) != 0L)
            stop("'strand' must be specified when 'rname' is not empty")
        strand <- strand()
    } else if (!is.factor(strand)
           || !identical(levels(strand), levels(strand())))
        stop("invalid 'strand' argument")
    if (!is.character(cigar) || any(is.na(cigar)))
        stop("'cigar' must be a character vector with no NAs")
    rglist <- cigarToIRangesListByAlignment(cigar, pos)
    strand <- .logicalAsCompactRawVector(strand == "-")
    new("Alignments0", rname=rname, strand=strand,
                       cigar=cigar, rglist=rglist)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("GappedAlignments", "Alignments0",
    function(from)
        Alignments0(rname=rname(from), strand=strand(from),
                    pos=start(from), cigar=cigar(from))

)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "Alignments0",
    function(x, i, j, ... , drop=TRUE)
    {
        i <- callNextMethod()
        x@cigar <- x@cigar[i]
        x@rname <- x@rname[i]
        x@strand <- .subsetCompactRawVector(x@strand, i)
        x@rglist <- x@rglist[i]
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "updateCigarAndStart" method.
###
### Performs atomic update of the cigar/start information.
###

setMethod("updateCigarAndStart", "Alignments0",
    function(x, cigar=NULL, start=NULL)
    {
        if (is.null(cigar))
            cigar <- cigar(x)
        else if (!is.character(cigar) || length(cigar) != length(x))
            stop("when not NULL, 'cigar' must be a character vector ",
                 "of the same length as 'x'")
        if (is.null(start))
            start <- start(x)
        else if (!is.integer(start) || length(start) != length(x))
            stop("when not NULL, 'start' must be an integer vector ",
                 "of the same length as 'x'")
        rglist <- cigarToIRangesListByAlignment(cigar, start)
        ## Atomic update (until the 2 slots are updated, x@cigar and x@rglist
        ## will be temporarily out of sync):
        x@cigar <- cigar
        x@rglist <- rglist
        x
    }
)

