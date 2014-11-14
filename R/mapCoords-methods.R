### =========================================================================
### 'mapCoords' and 'pmapCoords' methods
### -------------------------------------------------------------------------
###

### Generics are in IRanges.

.listCumsum <- function(x) {
    x_unlisted <- unlist(x, use.names=FALSE)
    x_cumsum <- cumsum(as.numeric(x_unlisted))
    x_part <- PartitioningByWidth(elementLengths(x))
    x_cumsum - rep(x_cumsum[start(x_part)] - x_unlisted[start(x_part)],
                   width(x_part))
}

.listCumsumShifted <- function(x) {
    cs <- .listCumsum(x)
    shifted <- c(0L, head(cs, -1))
    shifted[start(PartitioningByWidth(elementLengths(x)))] <- 0L
    shifted
}

## This method differs from sort() in that negative strand
## elements are returned highest value to lowest.
.orderElementsByTranscription <- function(x) {
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
    ## order by strand 
    neg <- strand(gr)[pstart] == "-"
    ord <- S4Vectors:::mseq(ifelse(neg, pend, pstart),
                            ifelse(neg, pstart, pend))
    res <- relist(gr[ord], x)
    res@unlistData$unordered <- original[idx[ord]] 
    res
}

.mapCoords <- function(from, to, ..., ignore.strand, elt.hits, p=FALSE) {
    if (ignore.strand)
        strand(to) <- "*"

    ## sort elements of 'to' by chrom, position and strand
    to <- .orderElementsByTranscription(to)
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
    cumsums <- .listCumsumShifted(width(to))
    cumPosition <- shift(eltPosition, 1L + cumsums[sHits])

    toInd <- togroup(to)[sHits]
    if (elt.hits)
        mcols <- DataFrame(queryHits=qHits, subjectHits=toInd, 
                           eltHits=mcols(gr)$unordered[subjectHits(ol)])
    else mcols <- DataFrame(queryHits=qHits, subjectHits=toInd)

    GRanges(seqnames(gr)[sHits], cumPosition, strand = strand(gr[sHits]), mcols)
}

### mapCoords:

setMethod("mapCoords", c("GenomicRanges", "GRangesList"), 
    function(from, to, ..., ignore.strand=TRUE, elt.hits=FALSE) 
        .mapCoords(from, to, ..., ignore.strand=ignore.strand,
                   elt.hits=elt.hits, p=FALSE)
)

setMethod("mapCoords", c("GenomicRanges", "GenomicRanges"), 
    function(from, to, ..., ignore.strand=TRUE, elt.hits=FALSE)
        callGeneric(from, relist(to, PartitioningByEnd(seq_along(to))),
                    ..., ignore.strand=ignore.strand, elt.hits=elt.hits)
)

### pmapCoords:

setMethod("pmapCoords", c("GenomicRanges", "GRangesList"), 
    function(from, to, ..., ignore.strand=TRUE, elt.hits=FALSE) 
    {
        if (length(from) != length(to))
            stop("'from' and 'to' must have the same length")
        ## FIXME: should be implemented as pfindOverlaps()
        .mapCoords(from, to, ..., ignore.strand=ignore.strand,
                   elt.hits=elt.hits, p=TRUE)
    }
)
