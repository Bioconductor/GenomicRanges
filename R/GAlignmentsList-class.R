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

###   ranges(x)   - IRangesList object of the same length as 'x'.
###   granges(x)  - GRangesList object of the same length as 'x'.
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

### Not implemented:
###   introns(x)  - Extract the N gaps in a GRangesList object of the same
###                 length as 'x'.
###   grglist(x)  - CompressedList (of GRangesList objects) the same 
###                 length as 'x'.
###   rglist(x)   - CompressedIRangesList object of the same length as 'x'.
###   No coercion to RangedDataList
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setMethod("ranges", "GAlignmentsList",
    function(x, ...) 
    {
        unlistData <- IRanges(start=x@unlistData@start, 
                          width=width(x@unlistData), 
                          names=x@unlistData@NAMES) 
        new2("CompressedIRangesList",
             unlistData=unlistData, partitioning=x@partitioning,
             elementMetadata=x@elementMetadata, check=FALSE)
    }
)

setMethod("granges", "GAlignmentsList",
    function(x) 
    {
        unlistData <- .GappedAlignmentsToGRanges(seqnames(x@unlistData), 
            start(x@unlistData), width(x@unlistData), strand(x@unlistData), 
            seqinfo(x@unlistData), names(x@unlistData))
        new2("GRangesList", unlistData=unlistData,
            partitioning=x@partitioning, elementMetadata=x@elementMetadata)
    }
)

setAs("GAlignmentsList", "RangesList", function(from) ranges(from))

setAs("GAlignmentsList", "GRangesList", function(from) granges(from))

setMethod("as.data.frame", "GAlignmentsList",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        if (missing(row.names))
            row.names <- names(x@unlistData)
        if (is.null(names(x)))
            element <- rep.int(seq_len(length(x)), elementLengths(x))
        else
            element <- rep.int(names(x), elementLengths(x))
        data.frame(element=element,
                   as.data.frame(unlist(x, use.names=FALSE),
                                 row.names=row.names),
                   stringsAsFactors=FALSE)
    }
)

.GAlignmentsListAsCompressedIRangesList <- function(from)
{
    ans_ranges <- IRanges(start=from@unlistData@start,
                          width=width(from@unlistData),
                          names=from@unlistData@NAMES)
    ans_ranges@elementMetadata <- from@unlistData@elementMetadata
    new("CompressedIRangesList",
        unlistData=ans_ranges,
        partitioning=from@partitioning,
        elementMetadata=from@elementMetadata)
}

setAs("GAlignmentsList", "CompressedIRangesList",
    .GAlignmentsListAsCompressedIRangesList
)

setAs("GAlignmentsList", "IRangesList",
    .GAlignmentsListAsCompressedIRangesList
)

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

