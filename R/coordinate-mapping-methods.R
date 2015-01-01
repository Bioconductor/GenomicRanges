### =========================================================================
### [p]mapToGenome() and [p]mapToTranscript() methods
### -------------------------------------------------------------------------
###

### I. Mapping strategy

### genome -> transcript:
### The *ToTranscript methods use findOverlaps for mapping. mapToTranscript 
### method attempts to map all elements of 'x' to all elements of 'alignment'
### where the parallel method only maps the i-th element of 'x' with the
### i-th element of 'alignment'.
###
### - findOverlaps handles strand
### - shift coords manually when 'alignment' is GRangesList
###
### transcript -> genome:
### The *ToGenome methods use the C code beneath transcriptLocs2refLocs()
### for mapping. mapToGenome does not attempt to map all i-j combinations but
### instead uses a name match between 'x' and 'alignment' to determine
### mapping pairs. The parallel pmapToGenome only attempts to map the i-th 
### element of 'x' with the i-th element of 'alignment'.
###
### - check strand manually
### - tlocs2rlocs handles shifting of coords when 'alignment' is GRangesList


### II. Output format

### Parallel methods return an object the same shape as 'x'. Ranges with
### a strand mismatch or non-hit are returned as zero-width ranges.
###
### Non-parallel methods return an object that varies in length like a
### Hits object. The result only contains mapped records, strand mismatch
### and non-hits are not returned.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generics 
###

setGeneric("mapToGenome", signature=c("x", "alignment"),
    function(x, alignment, ...) 
        standardGeneric("mapToGenome")
)

setGeneric("pmapToGenome", signature=c("x", "alignment"),
    function(x, alignment, ...) 
        standardGeneric("pmapToGenome")
)

setGeneric("mapToTranscript", signature=c("x", "alignment"),
    function(x, alignment, ...) 
        standardGeneric("mapToTranscript")
)

setGeneric("pmapToTranscript", signature=c("x", "alignment"),
    function(x, alignment, ...) 
        standardGeneric("pmapToTranscript")
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers
###

### 'x' is a GRangesList
### Returns a GRangesList with sorted elements.
### The method differs from sort() in that negative strand elements are
### returned highest value to lowest.
.orderElementsByTranscription <- function(x, ignore.strand) {
    original <- unlist(sapply(elementLengths(x), function(xx) 1:xx), 
                       use.names=FALSE)
    ## order by position
    gr <- unlist(x, use.names = FALSE)
    idx <- order(togroup(x), start(gr))
    gr <- gr[idx]
    part <- PartitioningByWidth(x)
    ## handle zero-width ranges
    pstart <- start(part)[width(part) != 0L]
    pend <- end(part)[width(part) != 0L]

    if (ignore.strand) {
        ord <- S4Vectors:::mseq(pstart, pend)
    } else {
        neg <- strand(gr)[pstart] == "-"
        ord <- S4Vectors:::mseq(ifelse(neg, pend, pstart),
                                ifelse(neg, pstart, pend))
    }
    res <- relist(gr[ord], x)
    res@unlistData$unordered <- original[idx[ord]] 
    res
}

### 'x' is an IntegerList or NumericList
### Returns a numeric vector of cumulative sums within list elements.
.listCumsumShifted <- function(x) {
    cs <- unlist(cumsum(x), use.names=FALSE)
    shifted <- c(0L, head(cs, -1))
    shifted[start(PartitioningByWidth(elementLengths(x)))] <- 0L
    shifted
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToTranscript() and pmapToTranscript() methods
###

.mapToTranscript <- function(x, alignment, hits, ignore.strand) { 
    if (length(hits)) {
        alignmentHits <- subjectHits(hits)
        xHits <- queryHits(hits)
        xrange <- ranges(x)[xHits]
        bounds <- ranges(alignment)[alignmentHits]

        ## location wrt to start of individual list elements
        if (ignore.strand) {
            xrange <- shift(xrange, - start(bounds) + 1L)
        } else {
            neg <- as.vector(strand(alignment)[alignmentHits] == "-")
            negstart <- end(bounds)[neg] - end(xrange)[neg] + 1L
            xrange[neg] <- IRanges(negstart, width=width(xrange)[neg])
            xrange[!neg] <- shift(xrange, - start(bounds) + 1L)
        }
        GRanges(seqnames(x)[xHits], xrange, strand(x)[xHits],
                DataFrame(xHits, alignmentHits))
    } else {
        ans <- GRanges()
        mcols(ans) <- DataFrame(xHits=integer(), alignmentHits=integer())
        ans
    }
}

.pmapToTranscript <- function(x, alignment, hits, ignore.strand) { 
    starts <- rep(1L, length(x))
    ends <- rep(0L, length(x))
    map <- .mapToTranscript(x, alignment, hits, ignore.strand)
    if (length(map)) {
        xHits <- mcols(map)$xHits
        starts[xHits] <- start(map)
        ends[xHits] <- end(map)
    }
    ## zero-width seqname
    seqname <- as.character(seqnames(x))
    seqname[setdiff(seq_along(x), queryHits(hits))] <- "unmapped"
    GRanges(seqname, IRanges(starts, ends, names=names(x)), strand(x))
}

setMethod("mapToTranscript", c("GenomicRanges", "GenomicRanges"), 
    function(x, alignment, ignore.strand=TRUE, ...)
    {
        hits <- findOverlaps(x, alignment, type="within", 
                             ignore.strand=ignore.strand)
        .mapToTranscript(x, alignment, hits, ignore.strand)
    }
)

setMethod("mapToTranscript", c("GenomicRanges", "GRangesList"), 
    function(x, alignment, ignore.strand=TRUE, ...) 
    {
        strand <- strand(x)
        if (!ignore.strand)
            if (!all(elementLengths(runLength(strand(alignment))) == 1))
                stop(paste0("when ignore.strand=TRUE all inner list ",
                            "elements of 'alignments' must be the same strand"))

        ## order within list elements by strand
        alignment <- .orderElementsByTranscription(alignment, ignore.strand)
        hits <- findOverlaps(x, unlist(alignment, use.names=FALSE), 
                             type="within", ignore.strand=ignore.strand)
        map <- .mapToTranscript(x, unlist(alignment, use.names=FALSE), 
                                hits, ignore.strand)

        ## location wrt start of concatenated list elements
        if (length(map)) {
            shifted <- .listCumsumShifted(width(alignment))
            map <- shift(map, shifted[subjectHits(hits)])
            mcols(map)$alignmentHits <- togroup(alignment)[subjectHits(hits)]
        }
        map
    }
)

setMethod("pmapToTranscript", c("GenomicRanges", "GenomicRanges"), 
    function(x, alignment, ignore.strand=TRUE, ...) 
    {
        if (length(x) != length(alignment))
            stop("'x' and 'alignment' must have the same length")

        hits <- findOverlaps(x, alignment, type="within", 
                             ignore.strand=ignore.strand)
        ith_hits <- queryHits(hits) == subjectHits(hits)
        hits <- hits[ith_hits] 
        .pmapToTranscript(x, alignment, hits, ignore.strand)
    }
)

setMethod("pmapToTranscript", c("GenomicRanges", "GRangesList"), 
    function(x, alignment, ignore.strand=TRUE, ...) 
    {
        if (length(x) != length(alignment))
            stop("'x' and 'alignment' must have the same length")
        if (!ignore.strand)
            if (!all(elementLengths(runLength(strand(alignment))) == 1))
                stop(paste0("when ignore.strand=TRUE all inner list ",
                            "elements of 'alignments' must be the same strand"))

        ## order within list elements
        alignment <- .orderElementsByTranscription(alignment, ignore.strand)
        hits <- findOverlaps(x, unlist(alignment, use.names=FALSE), 
                             type="within", ignore.strand=ignore.strand)
        ## filter hits on outer list element index
        ith_hits <- queryHits(hits) == togroup(alignment)[subjectHits(hits)]
        hits <- hits[ith_hits] 
        ans <- .pmapToTranscript(x, unlist(alignment, use.names=FALSE), 
                                 hits, ignore.strand)

        ## location wrt start of concatenated list elements
        ## shift zero-width ranges only
        if (any(hasWidth <- width(ans) != 0L)) {
            shifted <- .listCumsumShifted(width(alignment))
            ans[hasWidth] <- shift(ans[hasWidth], shifted[subjectHits(hits)])
        }

        ans 
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToGenome() and pmapToGenome() methods
###

.mapToGenome <- function(x, alignment, hits, ignore.strand)
{
    xHits <- queryHits(hits)
    alignmentHits <- subjectHits(hits)
    ## strand mismatch
    if (!ignore.strand) {
        strand_x <- strand(x)[xHits]
        strandAlign <- unlist(runValue(strand(alignment)), 
                              use.names=FALSE)[alignmentHits]
        strandMismatch <- 
            as.logical(strand_x == "+" & strandAlign == "-" |
                       strand_x == "-" & strandAlign == "+")
        strand <- as.character(strandAlign)[alignmentHits]
    } else {
        strandMismatch <- logical(length(xHits))
        strand <- rep("+", length(xHits))
    }
    xStart <- as.list(start(x)[xHits])
    xEnd <- as.list(end(x)[xHits])
    alignStart <- as.list(start(alignment)[alignmentHits])
    alignEnd <- as.list(end(alignment)[alignmentHits])

    starts <- unlist(.Call("tlocs2rlocs", xStart, alignStart, 
                           alignEnd, strand, FALSE, FALSE),
                     use.names=FALSE) 
    ends <- unlist(.Call("tlocs2rlocs", xEnd, alignStart, 
                         alignEnd, strand, FALSE, FALSE), 
                   use.names=FALSE) 

    seqname <- as.character(seqnames(x)[xHits])
    if (any(skip <- is.na(starts) | is.na(ends) | strandMismatch)) {
        starts[skip] <- 1L 
        ends[skip] <- 0L
        seqname[skip] <- "unmapped"
    }

    GRanges(Rle(seqname), IRanges(starts, ends, names=names(x)[xHits]),
            strand=strand)
}

setMethod("pmapToGenome", c("GenomicRanges", "GRangesList"), 
    function(x, alignment, ignore.strand=TRUE, ...) 
    {
        if (length(x) && length(alignment)) {
            if (length(x) != length(alignment))
                stop("'x' and 'alignment' must have the same length")
            if (!ignore.strand)
                if (!all(elementLengths(runLength(strand(alignment))) == 1))
                    stop(paste0("when ignore.strand=TRUE all inner list ",
                                "elements of 'alignments' must have the ",
                                "same strand"))

            ## order within list elements
            alignment <- .orderElementsByTranscription(alignment, ignore.strand)
            hits <- Hits(seq_along(x), seq_along(x), length(x), length(x))
            .mapToGenome(x, alignment, hits, ignore.strand)
        } else GRanges()
    }
)

setMethod("pmapToGenome", c("GenomicRanges", "GenomicRanges"),
    function(x, alignment, ignore.strand=TRUE, ...)
        pmapToGenome(x, splitAsList(alignment, seq_along(alignment)), 
                     ignore.strand=ignore.strand, ...)
)

setMethod("pmapToGenome", c("Ranges", "GenomicRanges"),
    function(x, alignment, ...) { 
        gr <- GRanges(seqnames(alignment), x)
        hits <- Hits(seq_along(x), seq_along(x), length(x), length(x))
        .mapToGenome(gr, alignment, hits, TRUE)
    }
)

setMethod("mapToGenome", c("GenomicRanges", "GRangesList"), 
    function(x, alignment, ignore.strand=TRUE, ...) 
    {
        if (length(x) && length(alignment)) {
            if (is.null(xNames <- names(x)) || 
                is.null(alignmentNames <- names(alignment)))
                stop ("both 'x' and 'alignment' must have names")
            if (!ignore.strand)
                if (!all(elementLengths(runLength(strand(alignment))) == 1))
                    stop(paste0("when ignore.strand=TRUE all inner list ",
                                "elements of 'alignments' must have the ",
                                "same strand"))

            match0 <- match(alignmentNames, alignmentNames)
            match1 <- match(xNames, alignmentNames)
            group0 <- splitAsList(seq_along(alignmentNames), match0)
            group1 <- group0[match(na.omit(match1), names(group0))]
            xHits <- rep(which(!is.na(match1)), elementLengths(group1))
            alignmentHits <- unlist(group1, use.names=FALSE)
            if (!length(xHits <- na.omit(xHits)))
                stop ("none of 'names(x)' are in 'names(alignment)'")

            ## order within list elements
            alignment <- .orderElementsByTranscription(alignment, ignore.strand)
            hits <- Hits(xHits, alignmentHits, length(x), length(alignment))
            ans <- .mapToGenome(x, alignment, hits, ignore.strand) 
            ## remove zero-width ranges, add mcols 
            df <- DataFrame(xHits, alignmentHits)
            mcols(ans) <- df
            if (any(hasWidth <- width(ans) != 0L))
                ans <- ans[hasWidth]
        } else {
            ans <- GRanges()
            mcols(ans) <- DataFrame(xHits=integer(), alignmentHits=integer())
        }

        ans
    }
)

setMethod("mapToGenome", c("GenomicRanges", "GenomicRanges"),
    function(x, alignment, ignore.strand=TRUE, ...) 
    {
        grl <- splitAsList(alignment, seq_along(alignment))
        names(grl) <- names(alignment)
        mapToGenome(x, grl, ignore.strand=ignore.strand, ...)
    }
)

setMethod("mapToGenome", c("Ranges", "GenomicRanges"),
    function(x, alignment, ...) { 
        gr <- GRanges(seqnames(alignment), x)
        ranges(mapToGenome(gr, alignment, TRUE, ...))
    }
)
