### =========================================================================
### 'mapCoords' and 'pmapCoords' methods
### -------------------------------------------------------------------------
###

### Generics are in IRanges.

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers
###

### 'x' is a GRangesList
### Returns a GRangesList with sorted elements. This method differs from 
### sort() in that "-" strand elements are returned highest value to lowest.
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

.mapCoords <- function(from, to, ..., ignore.strand, elt.hits, p=FALSE) {
    if (ignore.strand)
        strand(to) <- "*"

    ## sort elements of 'to' by chrom, position and strand
    to <- .orderElementsByTranscription(to, ignore.strand=ignore.strand)
    gr <- unlist(to, use.names = FALSE)

    ## overlaps
    ol <- findOverlaps(from, gr, type="within", ignore.strand=ignore.strand)
    if (p) {
        ith_hits <- queryHits(ol) == togroup(to)[subjectHits(ol)]
        ol <- ol[ith_hits] 
    }

    sHits <- subjectHits(ol)
    qHits <- queryHits(ol)
    eltPosition <- ranges(from)[qHits]
    bounds <- ranges(gr)[sHits]

    ## location wrt start of individual list elements
    if (ignore.strand) {
      eltPosition <- shift(eltPosition, - start(bounds))
    } else {
      neg <- as.vector(strand(gr)[sHits] == "-")
      eltPosition[!neg] <- shift(eltPosition[!neg], - start(bounds)[!neg])
      eltPosition[neg] <- IRanges(end(bounds)[neg] - end(eltPosition)[neg],
                                  width=width(eltPosition)[neg])
    }
    ## location wrt start of combined list elements (e.g., transcript-level)
    shifted <- .listCumsumShifted(width(to))
    cumPosition <- shift(eltPosition, 1L + shifted[sHits])

    toInd <- togroup(to)[sHits]
    if (elt.hits)
        mcols <- DataFrame(fromHits=qHits, toHits=toInd, 
                           eltHits=mcols(gr)$unordered[subjectHits(ol)])
    else mcols <- DataFrame(fromHits=qHits, toHits=toInd)

    GRanges(seqnames(gr)[sHits], cumPosition, strand = strand(gr[sHits]), mcols)
}

### mapCoords:

setMethod("mapCoords", c("GenomicRanges", "GRangesList"), 
    function(from, to, ..., ignore.strand=TRUE, elt.hits=FALSE) 
    {
        .Deprecated("mapToTranscripts", old="mapCoords")
        .mapCoords(from, to, ..., ignore.strand=ignore.strand,
                   elt.hits=elt.hits, p=FALSE)
    }
)

setMethod("mapCoords", c("GenomicRanges", "GenomicRanges"), 
    function(from, to, ..., ignore.strand=TRUE, elt.hits=FALSE) 
    {
        callGeneric(from, relist(to, PartitioningByEnd(seq_along(to))),
                    ..., ignore.strand=ignore.strand, elt.hits=elt.hits)
    }
)

### pmapCoords:

setMethod("pmapCoords", c("GenomicRanges", "GRangesList"), 
    function(from, to, ..., ignore.strand=TRUE, elt.hits=FALSE) 
    {
        .Deprecated("mapToTranscripts", old="mapCoords")
        if (length(from) != length(to))
            stop("'from' and 'to' must have the same length")
        ## FIXME: should be implemented as pfindOverlaps()
        .mapCoords(from, to, ..., ignore.strand=ignore.strand,
                   elt.hits=elt.hits, p=TRUE)
    }
)
