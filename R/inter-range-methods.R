### =========================================================================
### Inter-range methods
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### range()
###

.interIntervalSplitVariable <- function(x, ignore.strand=FALSE)
{
    f1 <- seqnames(x)
    runValue(f1) <- 3L * (as.integer(runValue(f1)) - 1L)
    if (ignore.strand)
       return(f1)
    f2 <- strand(x)
    runValue(f2) <- as.integer(runValue(f2)) - 1L
    f1 + f2
}

.interIntervalGenomicRanges2 <- function(x, FUN, ignore.strand=FALSE, ...)
{
    x <- clone(x)
    mcols(x) <- NULL
    f <- .interIntervalSplitVariable(x, ignore.strand=ignore.strand)
    xIRangesList <- split(unname(ranges(x)), f)
    ansIRangesList <- FUN(xIRangesList, ...)
    k <- elementLengths(ansIRangesList)
    x_seqlevels <- seqlevels(x)
    strand_levels <- levels(strand())
    ansIRangesList_names <- as.integer(names(ansIRangesList))
    i1 <- ansIRangesList_names %/% 3L + 1L
    ans_seqnames <- Rle(factor(x_seqlevels[i1], levels=x_seqlevels), k)
    if (ignore.strand) {
        ans_strand <- Rle(strand("*"), sum(k))
    } else {
        i2 <- ansIRangesList_names %% 3L + 1L
        ans_strand <- Rle(factor(strand_levels[i2], levels=strand_levels), k)
    }
    update(x, seqnames=ans_seqnames,
           ranges=unlist(ansIRangesList, use.names=FALSE),
           strand=ans_strand,
           elementMetadata=new("DataFrame", nrows=length(ans_seqnames)))
}

setMethod("range", "GenomicRanges",
    function(x, ..., ignore.strand=FALSE, na.rm=FALSE)
        .interIntervalGenomicRanges2(unname(c(x, ...)), range, 
                                     ignore.strand=ignore.strand)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reduce()
###

setMethod("reduce", "GenomicRanges",
    function(x, drop.empty.ranges=FALSE, min.gapwidth=1L,
             with.inframe.attrib=FALSE, ignore.strand=FALSE)
    {
        if (!identical(with.inframe.attrib, FALSE))
            stop("'with.inframe.attrib' argument not supported ",
                 "when reducing a GenomicRanges object")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand)
            strand(x) <- "*"
        .interIntervalGenomicRanges2(x, reduce,
                                     drop.empty.ranges=drop.empty.ranges,
                                     min.gapwidth=min.gapwidth)
    }
)

.interIntervalGenomicRanges <- function(x, FUN, ignore.strand=FALSE, ...)
{
    x <- clone(x)
    mcols(x) <- NULL
    if (ignore.strand)
       f <- paste(seqnames(x), Rle(factor("*"), length(x)), sep="\r")
    else
       f <- paste(seqnames(x), strand(x), sep="\r")
    xIRangesList <- split(unname(ranges(x)), f)
    ansIRangesList <- FUN(xIRangesList, ...)
    k <- elementLengths(ansIRangesList)
    splitListNames <- strsplit(names(ansIRangesList), split="\r", fixed=TRUE)
    listNameMatrix <- matrix(as.character(unlist(splitListNames)), nrow=2L)
    ansSeqnames <- Rle(factor(listNameMatrix[1L, ], levels=seqlevels(x)), k)
    ansStrand <- Rle(strand(listNameMatrix[2L, ]), k)
    update(x, seqnames=ansSeqnames,
           ranges=unlist(ansIRangesList, use.names=FALSE),
           strand=ansStrand,
           elementMetadata=new("DataFrame", nrows=length(ansSeqnames)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### gaps()
###

setMethod("gaps", "GenomicRanges",
    function(x, start=1L, end=seqlengths(x))
    {
        seqlevels <- seqlevels(x)
        if (!is.null(names(start)))
            start <- start[seqlevels]
        if (!is.null(names(end)))
            end <- end[seqlevels]
        start <- IRanges:::recycleVector(start, length(seqlevels))
        start <- rep(start, each=3L)
        end <- IRanges:::recycleVector(end, length(seqlevels))
        end <- rep(end, each=3L)
        .interIntervalGenomicRanges(x, gaps, start=start, end=end)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### disjoin()
###

setMethod("disjoin", "GenomicRanges",
    function(x, ignore.strand=FALSE)
        .interIntervalGenomicRanges2(x, disjoin, ignore.strand=ignore.strand)
)

applyOnRangesBySpace <- function(x, FUN, ..., ignore.strand = FALSE) {
  if (ignore.strand)
    f <- seqnames(x)
  else
    f <- paste(seqnames(x), strand(x), sep = "\r")
  xIRangesList <- split(unname(ranges(x)), f)
  ans <- FUN(xIRangesList, ...)
  if (is(ans, "List")) # values per range, otherwise assumed to be per-space
    ans <- unsplit(ans, f)
  ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isDisjoint()
###

setMethod("isDisjoint", "GenomicRanges",
    function(x, ignore.strand=FALSE)
    {
        all(applyOnRangesBySpace(x, isDisjoint, ignore.strand = ignore.strand))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### disjointBins()
###

setMethod("disjointBins", "GenomicRanges",
          function(x, ignore.strand = FALSE) {
            applyOnRangesBySpace(x, disjointBins,
                                 ignore.strand = ignore.strand)
          })

