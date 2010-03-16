### =========================================================================
### Alignments1 objects
### -------------------------------------------------------------------------
###

### Second GappedAlignments implementation: Alignments1
### The genomic ranges of the alignments are stored in the 'grglist' slot of
### type GRangesList.
setClass("Alignments1",
    contains="GappedAlignments",
    representation(
        grglist="GRangesList"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("rname", "Alignments1",
    function(x)
    {
        xgrg <- x@grglist
        as.factor(seqnames(xgrg@unlistData))[xgrg@partitioning@end]
    }
)

setReplaceMethod("rname", "Alignments1",
    function(x, value)
    {
        value <- normargRNameReplaceValue(x, value, ans.type="Rle")
        value <- rep.int(value, elementLengths(x@grglist))
        seqnames(x@grglist@unlistData) <- value
        x
    }
)

setMethod("strand", "Alignments1",
    function(x)
    {
        xgrg <- x@grglist
        as.factor(strand(xgrg@unlistData))[xgrg@partitioning@end]
    }
)

setMethod("grglist", "Alignments1", function(x) x@grglist)

setMethod("rglist", "Alignments1",
    function(x) as(ranges(x@grglist), "CompressedNormalIRangesList")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

Alignments1 <- function(rname=factor(), strand=NULL,
                        pos=integer(), cigar=character())
{
    rglist <- cigarToIRangesListByAlignment(cigar, pos)
    if (is.null(strand)) {
        if (length(rname) != 0L)
            stop("'strand' must be specified when 'rname' is not empty")
        strand <- strand()
    } else if (!is.factor(strand)
            || !identical(levels(strand), levels(strand())))
        stop("invalid 'strand' argument")
    grglist <- GappedAlignmentsAsGRangesList(rname, strand, rglist)
    new("Alignments1", cigar=cigar, grglist=grglist)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("GappedAlignments", "Alignments1",
    function(from)
        Alignments1(rname=rname(from), strand=strand(from),
                    pos=start(from), cigar=cigar(from))

)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "Alignments1",
    function(x, i, j, ... , drop=TRUE)
    {
        i <- callNextMethod()
        x@cigar <- x@cigar[i]
        x@grglist <- x@grglist[i]
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "updateCigarAndStart" method.
###
### Performs atomic update of the cigar/start information.
###

setMethod("updateCigarAndStart", "Alignments1",
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
        grglist <- GappedAlignmentsAsGRangesList(rname(x), strand(x), rglist)
        ## Atomic update (until the 2 slots are updated, x@cigar and x@grglist
        ## will be temporarily out of sync):
        x@cigar <- cigar
        x@grglist <- grglist
        x
    }
)

