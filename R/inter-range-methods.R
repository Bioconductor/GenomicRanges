### =========================================================================
### Inter-range methods
### -------------------------------------------------------------------------
###

### The methods documented in this page on this page are consistent with
### those in IRanges inter-range-methods.R
### range()
### reduce()
### gaps()
### disjoin()
### isDisjoint()
### disjointBins()


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

### 'revmap_unlisted' and 'revmap_partitioning': IntegerList and Partitioning
### objects representing the "revmap object", which is *conceptually* an
### IntegerListList object (of the same length as 'revmap_partitioning').
### *Conceptually* because, well, we don't actually have such container...
### 'old2new': IntegerList of the same length as the "revmap object".
.translateRevmap <- function(revmap_unlisted, revmap_partitioning, old2new)
{
    ## 'times' has the length of the "revmap object".
    times <- sum(relist(width(PartitioningByEnd(revmap_unlisted)),
                        revmap_partitioning))
    ## 'offset' has the length of 'revmap_unlisted@unlistData'.
    offset <- rep.int(start(PartitioningByEnd(old2new)) - 1L, times)
    revmap_flat <- revmap_unlisted@unlistData
    revmap_unlisted@unlistData <- old2new@unlistData[revmap_flat + offset]
    revmap_unlisted
}

.interIntervalGenomicRanges2 <- function(x, FUN,
                                         with.revmap=FALSE, with.mapping=FALSE,
                                         ignore.strand=FALSE, ...)
{
    x <- clone(x)
    mcols(x) <- NULL
    f <- .interIntervalSplitVariable(x, ignore.strand=ignore.strand)
    xIRangesList <- split(unname(ranges(x)), f)
    if (with.revmap || with.mapping) {
        ansIRangesList <- FUN(xIRangesList,
                              with.revmap=with.revmap,
                              with.mapping=with.mapping, ...)
    } else {
        ansIRangesList <- FUN(xIRangesList, ...)
    }
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
    ans_ranges <- unlist(ansIRangesList, use.names=FALSE)
    ans_mcols <- mcols(ans_ranges)
    if (is.null(ans_mcols)) {
        ans_mcols <- new("DataFrame", nrows=length(ans_seqnames))
    } else {
        mcols(ans_ranges) <- NULL
        if (with.revmap || with.mapping) {
            if (with.revmap) {
                revmap_colname <- "revmap"
            } else {
                revmap_colname <- "mapping"
            }
            revmap_unlisted <- ans_mcols[[revmap_colname]]
            revmap_partitioning <- PartitioningByEnd(ansIRangesList)
            old2new <- splitAsList(seq_along(x), f)
            revmap_unlisted2 <- .translateRevmap(revmap_unlisted,
                                                 revmap_partitioning,
                                                 old2new)
            ans_mcols[[revmap_colname]] <- revmap_unlisted2
        }
    }
    if (length(extraColumnSlotNames(x)) > 0L)
        x <- new("GRanges", x)
    update(x, seqnames=ans_seqnames,
              ranges=ans_ranges,
              strand=ans_strand,
              elementMetadata=ans_mcols)
}

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
    if (length(extraColumnSlotNames(x)) > 0L)
        x <- new("GRanges", x)
    update(x, seqnames=ansSeqnames,
           ranges=unlist(ansIRangesList, use.names=FALSE),
           strand=ansStrand,
           elementMetadata=new("DataFrame", nrows=length(ansSeqnames)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### range()
###

setMethod("range", "GenomicRanges",
    function(x, ..., ignore.strand=FALSE, na.rm=FALSE)
        .interIntervalGenomicRanges2(unname(c(x, ...)), range, 
                                     ignore.strand=ignore.strand)
)

### TODO: Support the 'ignore.strand' arg.
setMethod("range", "GRangesList",
    function(x, ..., na.rm=FALSE)
    {
        if (length(list(...)) != 0L)
            stop("\"range\" method for GRangesList objects only ",
                 "takes a single object")
        gr <- deconstructGRLintoGR(x)
        ## "range" method for GRanges objects is fast.
        gr <- range(gr)
        reconstructGRLfromGR(gr, x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### reduce()
###

setMethod("reduce", "GenomicRanges",
    function(x, drop.empty.ranges=FALSE, min.gapwidth=1L,
                with.revmap=FALSE, with.mapping=FALSE,
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
                                     min.gapwidth=min.gapwidth,
                                     with.revmap=with.revmap,
                                     with.mapping=with.mapping)
    }
)

### TODO: Support the 'ignore.strand' arg.
setMethod("reduce", "GRangesList",
    function(x, drop.empty.ranges=FALSE, min.gapwidth=1L,
             with.inframe.attrib=FALSE)
    {
        if (!identical(with.inframe.attrib, FALSE)) 
            stop("'with.inframe.attrib' argument is not supported ", 
                 "when reducing a GRangesList object")
        gr <- deconstructGRLintoGR(x)
        ## "reduce" method for GRanges objects is fast.
        gr <- reduce(gr, drop.empty.ranges=drop.empty.ranges,
                         min.gapwidth=min.gapwidth)
        reconstructGRLfromGR(gr, x)
    }
)


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
        start <- S4Vectors:::recycleVector(start, length(seqlevels))
        start <- rep(start, each=3L)
        end <- S4Vectors:::recycleVector(end, length(seqlevels))
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

setMethod("disjoin", "GRangesList",
    function(x, ...)
{
    gr <- deconstructGRLintoGR(x)
    d <- disjoin(gr, ...)
    reconstructGRLfromGR(d, x)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isDisjoint()
###

setMethod("isDisjoint", "GenomicRanges",
    function(x, ignore.strand=FALSE)
    {
        all(applyOnRangesBySpace(x, isDisjoint, ignore.strand = ignore.strand))
    }
)

setMethod("isDisjoint", "GRangesList",
    function(x, ignore.strand = FALSE)
    {
        gr <- deconstructGRLintoGR(x)

        if (ignore.strand) 
            xIRangesList <- split(unname(ranges(gr)), paste(seqnames(gr),
                           Rle(factor(rep("+", length(gr)))), sep = "\r"))
        else 
            xIRangesList <- split(unname(ranges(gr)),
                                  paste(seqnames(gr), strand(gr), sep = "\r"))
 
        ansIRanges <- isDisjoint(xIRangesList)
        splitListNames <- strsplit(names(ansIRanges), split="\r")
        snames <- strsplit(unlist(lapply(splitListNames, "[[", 1L)), "|",
                           fixed=TRUE)
        m12 <- matrix(as.integer(unlist(snames)), ncol=2, byrow=TRUE)

        ansIRangesList <- split(ansIRanges,
                                factor(m12[, 1L], levels=seq_len(length(x))))
        ans <-  unlist(lapply(ansIRangesList, all))
        names(ans) <- names(x)
        ans
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

