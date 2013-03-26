### =========================================================================
### GAlignmentsList objects
### -------------------------------------------------------------------------
###

setClass("GAlignmentsList",
    contains="CompressedList",
    representation(
        unlistData="GappedAlignments",
        elementMetadata="DataFrame"
    ),
    prototype(
        elementType="GappedAlignments"
    )
)

### Formal API:
###   names(x)    - NULL or character vector.
###   length(x)   - single integer. Nb of alignments in 'x'.
###   seqnames(x) - 'factor' Rle of the same length as 'x'.
###   rname(x)    - same as 'seqnames(x)'.
###   seqnames(x) <- value - replacement form of 'seqnames(x)'.
###   rname(x) <- value - same as 'seqnames(x) <- value'.
###   cigar(x)    - character vector of the same length as 'x'.
###   strand(x)   - 'factor' Rle of the same length as 'x' (levels: +, -, *).
###   qwidth(x)   - integer vector of the same length as 'x'.
###   start(x), end(x), width(x) - integer vectors of the same length as 'x'.
###   ngap(x)     - integer vector of the same length as 'x'.

###   grglist(x)  - GRangesList object of the same length as 'x'.
###   granges(x)  - GRanges object of the same length as 'x'.
###   introns(x)  - Extract the N gaps in a GAlignmentsList object of the same
###                 length as 'x'.
###   rglist(x)   - CompressedIRangesList object of the same length as 'x'.
###   ranges(x)   - IRanges object of the same length as 'x'.
###   as.data.frame(x) - data.frame with 1 row per alignment in 'x'.

###   show(x)     - compact display in a data.frame-like fashion.
###   GAlignmentsList(x, ...) - constructor.
###   x[i]        - GAlignmentsList object of the same class as 'x'
###                 (endomorphism).
###
###   findOverlaps(query, subject) - 'query' or 'subject' or both are
###                 GappedAlignments objects. Just a convenient wrapper for
###                 'findOverlaps(grglist(query), subject, ...)', etc...
###
###   countOverlaps(query, subject) - 'query' or 'subject' or both are
###                 GappedAlignments objects. Just a convenient wrapper for
###                 'countOverlaps(grglist(query), subject, ...)', etc...
###
###   subsetByOverlaps(query, subject) - 'query' or 'subject' or both are
###                 GappedAlignments objects.
###

###   qnarrow(x, start=NA, end=NA, width=NA) - GAlignmentsList object of the
###                 same length and class as 'x' (endomorphism).
###
###   narrow(x, start=NA, end=NA, width=NA) - GAlignmentsList object of the
###                 same length and class as 'x' (endomorphism).
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setMethod("seqnames", "GAlignmentsList", 
    function(x) 
        new2("CompressedRleList",
             unlistData=x@unlistData@seqnames, 
             partitioning=x@partitioning, check=FALSE)
)

setMethod("rname", "GAlignmentsList", 
    function(x) 
        new2("CompressedRleList",
             unlistData=x@unlistData@seqnames, 
             partitioning=x@partitioning, check=FALSE)
)

setMethod("cigar", "GAlignmentsList", 
    function(x) 
        new2("CompressedCharacterList",
             unlistData=x@unlistData@cigar, 
             partitioning=x@partitioning, check=FALSE)
)

setMethod("strand", "GAlignmentsList",
    function(x)
        new2("CompressedRleList",
             unlistData=x@unlistData@strand, 
             partitioning=x@partitioning, check=FALSE)
)

setMethod("qwidth", "GAlignmentsList",
    function(x)
        new2("CompressedIntegerList",
             unlistData=cigarToWidth(x@unlistData@cigar),
             partitioning=x@partitioning, check=FALSE)
)

setMethod("ngap", "GAlignmentsList",
    function(x)
        new2("CompressedIntegerList",
             unlistData=unname(elementLengths(rglist(x@unlistData))) - 1L,
             partitioning=x@partitioning, check=FALSE)
)

setMethod("start", "GAlignmentsList",
    function(x, ...)
        new2("CompressedIntegerList",
             unlistData=x@unlistData@start,
             partitioning=x@partitioning, check=FALSE)
)

setMethod("end", "GAlignmentsList",
    function(x, ...)
        new2("CompressedIntegerList",
             unlistData=end(x@unlistData),
             partitioning=x@partitioning, check=FALSE)
)

setMethod("width", "GAlignmentsList",
    function(x)
        new2("CompressedIntegerList",
             unlistData=cigarToWidth(x@unlistData@cigar),
             partitioning=x@partitioning, check=FALSE)
)

setMethod("seqinfo", "GAlignmentsList", function(x) seqinfo(x@unlistData))

setMethod("elementMetadata", "GAlignmentsList", .getElementMetadataList)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###

setReplaceMethod("names", "GAlignmentsList",
    function(x, value) 
    {
        if (!is(value, "AtomicList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an AtomicList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        if (!is(value, "Rle"))
            value <- Rle(factor(value))
        else if (!is.factor(runValue(value)))
            runValue(value) <- factor(runValue(value))
        names(x@unlistData) <- value
        x
    }
)

setReplaceMethod("rname", "GAlignmentsList",
    function(x, value) `seqnames<-`(x, value)
)

setReplaceMethod("elementMetadata", "GAlignmentsList", 
    .replaceElementMetadataList
)

setReplaceMethod("strand", "GAlignmentsList", .replaceStrandList)

setReplaceMethod("seqinfo", "GAlignmentsList", .replaceSeqinfoList)

setReplaceMethod("seqnames", "GAlignmentsList", .replaceSeqnamesList)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GAlignmentsList.mcols <- function(x)
{
}

.valid.GAlignmentsList <- function(x)
{
    c(.valid.GAlignmentsList.mcols(x))
}

setValidity2("GAlignmentsList", .valid.GAlignmentsList)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

GAlignmentsList <- function(...)
{
    listData <- list(...)
    if (length(listData) == 0L) {
        unlistData <- GappedAlignments()
    } else {
        if (length(listData) == 1L && is.list(listData[[1L]]))
            listData <- listData[[1L]]
        if (!all(sapply(listData, is, "GappedAlignments")))
            stop("all elements in '...' must be GappedAlignments objects")
        unlistData <- suppressWarnings(do.call("c", unname(listData)))
    }
    relist(unlistData, PartitioningByEnd(listData))
}

makeGAlignmentsListFromFeatureFragments <- function(seqnames=Rle(factor()),
                                                fragmentPos=list(),
                                                fragmentCigar=list(),
                                                strand=character(0),
                                                sep=",")
{
    fragmentPos <- .normargListOfIntegers(fragmentPos, sep, "fragmentPos")
    nfrag_per_feature <- elementLengths(fragmentPos)
    pos <- unlist(fragmentPos, recursive=FALSE, use.names=FALSE)

    ncigar_per_elt <- elementLengths(fragmentCigar)
    if (length(ncigar_per_elt) != 0L) {
        if (length(nfrag_per_feature) == 0L)
            nfrag_per_feature <- ncigar_per_elt
        else if (!identical(ncigar_per_elt, nfrag_per_feature))
            stop("'fragmentPos' and 'fragmentCigar' have ",
                 "incompatible \"shapes\"")
    }
    cigar <- unlist(fragmentCigar, recursive=FALSE, use.names=FALSE)

    nfrag <- sum(nfrag_per_feature)
    if (nfrag == 0L) {
        ## Cannot blindly subset by FALSE because it doesn't work on a
        ## zero-length Rle.
        if (length(seqnames) != 0L)
            seqnames <- seqnames[FALSE]
        if (length(strand) != 0L)
            strand <- strand[FALSE]
    } else {
        if (length(seqnames) != length(nfrag_per_feature) ||
            length(strand) != length(nfrag_per_feature))
            stop("length of 'seqnames' and/or 'strand' is incompatible ",
                 "with fragmentPos/Cigar")
        seqnames <- rep.int(seqnames, nfrag_per_feature)
        strand <- rep.int(strand, nfrag_per_feature)
    }
    unlistData <- GappedAlignments(seqnames=seqnames, pos=pos, cigar=cigar,
                                   strand=strand)
    partitioning <- PartitioningByEnd(cumsum(nfrag_per_feature), names=NULL)
    relist(unlistData, partitioning)
}

setMethod("updateObject", "GAlignmentsList",
    function(object, ..., verbose=FALSE)
    {
        if (verbose)
            message("updateObject(object = 'GAlignmentsList')")
        if (is(try(validObject(object@unlistData, complete=TRUE), silent=TRUE),
               "try-error")) {
            object@unlistData <- updateObject(object@unlistData)
            return(object)
        }
        object
    }
)

readGAlignmentsList <- function(file, format="BAM", use.names=FALSE, ...)
{
    if (!isSingleString(format))
        stop("'format' must be a single string")
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (format == "BAM") {
        suppressMessages(library("Rsamtools"))
        ans <- readBamGAlignmentsList(file=file, use.names=TRUE, ...) 
        return(ans)
    }
    stop("only BAM format is supported at the moment")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setMethod("grglist", "GAlignmentsList",
    function(x, order.as.in.query=FALSE, drop.D.ranges=FALSE,
             ignore.strand=FALSE) 
    {
        if (ignore.strand)
            strand(x@unlistData) <- "*"
        rgl <- rglist(x@unlistData,
                      order.as.in.query=order.as.in.query,
                      drop.D.ranges=drop.D.ranges)
        eltlen <- elementLengths(rgl)
        gr <- GRanges(rep.int(seqnames(x@unlistData), eltlen), 
                      ranges=unlist(rgl), 
                      strand=rep.int(strand(x@unlistData), eltlen),
                      seqinfo=seqinfo(x@unlistData))
        relist(gr, .shiftPartition(x@partitioning, rgl@partitioning))
    }
)
 
setMethod("granges", "GAlignmentsList",
    function(x, ignore.strand=FALSE) 
    {
        if (ignore.strand)
            strand(x@unlistData) <- "*"
        gr <- range(grglist(x))@unlistData
        if (length(gr) != length(x))
            warning("Some ranges could not be combined. ",
                    "Consider using 'ignore.strand=TRUE'.")
        gr
    }
)

setMethod("introns", "GAlignmentsList",
    function(x, ignore.strand=FALSE)
    {
        if (ignore.strand)
            strand(x@unlistData) <- "*"
        grl <- introns(x@unlistData)
        f <- rep(seq_along(x@partitioning), width(x@partitioning))
        relist(grl@unlistData, 
               .shiftPartition(x@partitioning, grl@partitioning))
    }
)

setMethod("rglist", "GAlignmentsList",
    function(x, order.as.in.query=FALSE, drop.D.ranges=FALSE)
    {
        rgl <- rglist(x@unlistData)
        partitioning <- .shiftPartition(x@partitioning, rgl@partitioning) 
        relist(rgl@unlistData, partitioning) 
    }
)

## Adjust the widths of 'partition1' to accomodate the 
## increased / decreased number of splits in 'partition2'.
## The return value is a PartitioningByEnd object the same 
## length as 'partition1'.
.shiftPartition <- function(partition1, partition2)
{
    f <- rep(seq_along(partition1), width(partition1))
    w <- sapply(split(width(partition2), f), sum)
    if (any(w < 0))
        stop("width of 'partition2' cannot be negative")
    wdiff <- w - width(partition1)
    PartitioningByEnd(end(partition1) + cumsum(wdiff), names=names(partition1))
}

setMethod("ranges", "GAlignmentsList",
    function(x) 
    {
        rgl <- rglist(x@unlistData)
        lst <- relist(unlist(rgl), 
                      .shiftPartition(x@partitioning, rgl@partitioning))
        unlist(range(lst), use.names=FALSE)
    }
)

setAs("GAlignmentsList", "GRangesList", 
    function(from) grglist(from, ignore.strand=ignore.strand)
)
setAs("GAlignmentsList", "GRanges", 
    function(from) granges(from, ignore.strand=ignore.strand)
)
setAs("GAlignmentsList", "RangesList", 
    function(from) rglist(from, ignore.strand=ignore.strand)
)
setAs("GAlignmentsList", "Ranges", 
    function(from) ranges(from, ignore.strand=ignore.strand)
)
setMethod("as.data.frame", "GAlignmentsList", .GRangesListAsdataframe)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

setMethod("[", "GAlignmentsList", .sBracketSubsetGRList)
setReplaceMethod("[", "GAlignmentsList", .sBracketReplaceGRList)
setReplaceMethod("[[", "GAlignmentsList", .dBracketReplaceGRList)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

setMethod("show", "GAlignmentsList",
    function(object)
        .showList(object, showGappedAlignments, FALSE)
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "qnarrow" and "narrow" methods.
###

setMethod("qnarrow", "GAlignmentsList",
    function(x, start=NA, end=NA, width=NA) 
    {
        gal <- qnarrow(x@unlistData, start=start, end=end, width=width)
        relist(gal, x@partitioning)
    }
)

setMethod("narrow", "GAlignmentsList",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        gal <- narrow(x@unlistData, start=start, end=end, width=width,
                      use.names=use.names)
        relist(gal, x@partitioning)
    }
)
