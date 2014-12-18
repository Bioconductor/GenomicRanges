### =========================================================================
### 'mapToGenome' and 'mapToTranscript' methods
### -------------------------------------------------------------------------
###

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

### This method differs from sort() in that negative strand
### elements are returned highest value to lowest.
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
        ## order by strand 
        neg <- strand(gr)[pstart] == "-"
        ord <- S4Vectors:::mseq(ifelse(neg, pend, pstart),
                                ifelse(neg, pstart, pend))
    }
    res <- relist(gr[ord], x)
    res@unlistData$unordered <- original[idx[ord]] 
    res
}

### 'x' is an IntegerList or NumericList
### returns a numeric vector of cumulative sums within list elements
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
        position <- ranges(x)[xHits]
        bounds <- ranges(alignment)[alignmentHits]

        ## location wrt to start of individual list elements
        if (ignore.strand) {
            position <- shift(position, - start(bounds))
        } else {
            neg <- as.vector(strand(alignment)[alignmentHits] == "-")
            position[!neg] <- shift(position[!neg], - start(bounds)[!neg])
            position[neg] <- IRanges(end(bounds)[neg] - end(position)[neg],
                                        width=width(position)[neg])
        }
        ans <- GRanges(seqnames(x)[xHits], position, strand(x)[xHits],
                       DataFrame(xHits, alignmentHits))
    } else {
        ans <- GRanges()
        mcols(ans) <- DataFrame(xHits=integer(), alignmentHits=integer())
    }
    ans
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
        map <- .mapToTranscript(x, unlist(alignment, use.names=FALSE), hits,
                                ignore.strand)

        ## location wrt start of concatenated list elements
        if (length(map)) {
            shifted <- .listCumsumShifted(width(alignment))
            map <- shift(map, 1L + shifted[subjectHits(hits)])
            mcols(map)$alignmentHits <- togroup(alignment)[subjectHits(hits)]
        }
        map
    }
)

.pmapToTranscript <- function(x, alignment, hits, ignore.strand=TRUE) { 
    ## FIXME: handle non-hits re: pintersect arg 'resolve.empty'
    starts <- rep(1L, length(x))
    ends <- rep(0L, length(x))
    map <- .mapToTranscript(x, alignment, hits, ignore.strand)
    if (length(map)) {
        xHits <- mcols(map)$xHits
        starts[xHits] <- start(map)
        ends[xHits] <- end(map)
    }
    ans <- GRanges(seqnames(x), IRanges(starts, ends), strand(x))
    names(ans) <- names(x)
    ans
}

setMethod("pmapToTranscript", c("GenomicRanges", "GenomicRanges"), 
    function(x, alignment, ignore.strand=TRUE, ...) 
    {
        if (length(x) != length(alignment))
            stop("'x' and 'alignment' must have the same length")

        ## FIXME: pfindOverlaps?
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
        strand <- strand(x)
        if (!ignore.strand)
            if (!all(elementLengths(runLength(strand(alignment))) == 1))
                stop(paste0("when ignore.strand=TRUE all inner list ",
                            "elements of 'alignments' must be the same strand"))

        ## order within list elements by strand
        alignment <- .orderElementsByTranscription(alignment, ignore.strand)
        ## FIXME: pfindOverlaps?
        hits <- findOverlaps(x, unlist(alignment, use.names=FALSE), 
                             type="within", ignore.strand=ignore.strand)
        ith_hits <- queryHits(hits) == togroup(alignment)[subjectHits(hits)]
        hits <- hits[ith_hits] 
        map <- .pmapToTranscript(x, unlist(alignment, use.names=FALSE), 
                                 hits, ignore.strand)

        ## location wrt start of concatenated list elements
        ## shift zero-width ranges only
        if (any(hasWidth <- width(map) != 0L)) {
            shifted <- .listCumsumShifted(width(alignment))
            map[hasWidth] <- shift(map[hasWidth], 
                                   1L + shifted[subjectHits(hits)])
        }
        map
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToGenome() and pmapToGenome() methods
###


