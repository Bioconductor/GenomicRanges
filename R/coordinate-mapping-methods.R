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

## This method differs from sort() in that negative strand
## elements are returned highest value to lowest.
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

.listCumsumShifted <- function(x) {
  cs <- unlist(cumsum(x), use.names=FALSE)
  shifted <- c(0L, head(cs, -1))
  shifted[start(PartitioningByWidth(elementLengths(x)))] <- 0L
  shifted
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToTranscript() and pmapToTranscript() methods
###

### 'alignment' is a GRanges
.mapToTranscript_GR <- function(x, alignment, ignore.strand=TRUE) {
    ol <- findOverlaps(x, alignment, type="within", ignore.strand=ignore.strand)
    if (length(ol)) {
        sHits <- subjectHits(ol)
        qHits <- queryHits(ol)
        eltPosition <- ranges(x)[qHits]
        bounds <- ranges(alignment)[sHits]

        ## location wrt start of individual list elements
        if (ignore.strand) {
          eltPosition <- shift(eltPosition, - start(bounds))
        } else {
          neg <- as.vector(strand(alignment)[sHits] == "-")
          eltPosition[!neg] <- shift(eltPosition[!neg], - start(bounds)[!neg])
          eltPosition[neg] <- IRanges(end(bounds)[neg] - end(eltPosition)[neg],
                                      width=width(eltPosition)[neg])
        }
        ans <- eltPosition 
        mcols(ans) <- DataFrame(xHits=qHits, alignmentHits=sHits)
    } else {
        ans <- IRanges()
        mcols(ans) <- DataFrame(xHits=integer(), alignmentHits=integer())
    }
    ans
}

setMethod("mapToTranscript", c("GenomicRanges", "GenomicRanges"), 
    function(x, alignment, ignore.strand=TRUE, ...)
    {
        map <- .mapToTranscript_GR(x, alignment, ignore.strand) 
        GRanges(seqnames(x)[mcols(map)$xHits], map,
                strand(alignment[mcols(map)$xHits]))
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
        map <- .mapToTranscript_GR(x, 
                                   unlist(alignment, use.names=FALSE), 
                                   ignore.strand)

        if (length(map)) {
            xHits <- mcols(map)$xHits
            alignmentHits <- mcols(map)$alignmentHits
            ## adjust range positions to start of combined list elements
            shifted <- .listCumsumShifted(width(alignment))
            cumPosition <- shift(map, 1L + shifted[alignmentHits])
            GRanges(seqnames(x)[xHits], cumPosition, strand[xHits],
                    DataFrame(xHits, 
                              alignmentHits=togroup(alignment)[alignmentHits])) 
        } else {
            ans <- GRanges()
            mcols(ans) <- DataFrame(xHits=integer(), alignmentHits=integer())
            ans
        }
    }
)

### 'alignment' is a GRanges
.pmapToTranscript_GR <- function(x, alignment, ignore.strand, elts) {
    ol <- findOverlaps(x, alignment, type="within", ignore.strand=ignore.strand)
    ith_hits <- queryHits(ol) == elts[subjectHits(ol)]
    ol <- ol[ith_hits] 

    starts <- rep(1L, length(x))
    ends <- rep(0L, length(x))
    df <- DataFrame(xHits=seq_along(x), alignmentHits=seq_along(x))
    if (length(ol)) {
        sHits <- subjectHits(ol)
        qHits <- queryHits(ol)
        eltPosition <- ranges(x)[qHits]
        bounds <- ranges(alignment)[sHits]

        ## location wrt start of individual list elements
        if (ignore.strand) {
          eltPosition <- shift(eltPosition, - start(bounds))
        } else {
          neg <- as.vector(strand(alignment)[sHits] == "-")
          eltPosition[!neg] <- shift(eltPosition[!neg], - start(bounds)[!neg])
          eltPosition[neg] <- IRanges(end(bounds)[neg] - end(eltPosition)[neg],
                                      width=width(eltPosition)[neg])
        }
        ## FIXME: handle non-hits re: pintersect arg 'resolve.empty'
        starts[qHits] <- start(eltPosition)
        ends[qHits] <- end(eltPosition)
        df$alignmentHits[qHits] <- sHits
    }
    ans <- IRanges(starts, ends)
    mcols(ans) <- df
    ans
}

setMethod("pmapToTranscript", c("GenomicRanges", "GenomicRanges"), 
    function(x, alignment, ignore.strand=TRUE, ...) 
    {
        if (length(x) != length(alignment))
            stop("'x' and 'alignment' must have the same length")
        map <- .pmapToTranscript_GR(x, alignment, ignore.strand,
                                    togroup(alignment))
        mcols(map) <- NULL
        GRanges(seqnames(x)[mcols(map)$xHits], map,
                strand(alignment[mcols(map)$xHits]))
    }
)

setMethod("pmapToTranscript", c("GenomicRanges", "GRangesList"), 
    function(x, alignment, ignore.strand=TRUE, ...) 
    {
        strand <- strand(x)
        if (!ignore.strand)
            if (!all(elementLengths(runLength(strand(alignment))) == 1))
                stop(paste0("when ignore.strand=TRUE all inner list ",
                            "elements of 'alignments' must be the same strand")) 
        if (length(x) != length(alignment))
            stop("'x' and 'alignment' must have the same length")

        ## order within list elements by strand
        alignment <- .orderElementsByTranscription(alignment, ignore.strand)
        map <- .pmapToTranscript_GR(x, 
                                    unlist(alignment, use.names=FALSE),
                                    ignore.strand, togroup(alignment))

        ## adjust range positions to start of combined list elements
        ## shift only non-zero width ranges
        if (any(hasWidth <- width(map) != 0L)) {
            shifted <- .listCumsumShifted(width(alignment))
            alignmentHits <- mcols(map)$alignmentHits
            map[hasWidth] <- shift(map[hasWidth], 
                                   1L + shifted[alignmentHits[hasWidth]])
        }
        mcols(map) <- NULL
        GRanges(seqnames(x), map, strand(x))
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### mapToGenome() and pmapToGenome() methods
###
