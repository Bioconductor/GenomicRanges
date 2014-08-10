### =========================================================================
### mapCoords methods
### -------------------------------------------------------------------------
###

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

.orderElementsByTranscription <- function(x) {
  gr <- unlist(x, use.names = FALSE)
  gr <- gr[order(togroup(x), start(gr))]
  part <- PartitioningByWidth(x)
  ## handle zero-width ranges
  pstart <- start(part)[width(part) != 0L]
  pend <- end(part)[width(part) != 0L]
  neg <- strand(gr)[pstart] == "-"
  ord <- S4Vectors:::mseq(ifelse(neg, pend, pstart),
                          ifelse(neg, pstart, pend))
  relist(gr[ord], x)
}

setMethod("map", c("GenomicRanges", "GRangesList"), function(from, to) {
  .Defunct(msg="map() is defunct. Use mapCoords() instead.")
})

setMethod("mapCoords", c("GenomicRanges", "GRangesList"), 
  function(x, to, ..., ignore.strand = FALSE, elt.loc = FALSE,
           elt.hits = FALSE) {
    if (ignore.strand)
        strand(to@unlistData) <- "*"
    ## make sure 'to' is properly sorted by strand
    to <- .orderElementsByTranscription(to)

    ## find overlaps
    gr <- unlist(to, use.names = FALSE)
    ol <- findOverlaps(x, gr, type = "within", ..., 
                       ignore.strand=ignore.strand)
    shits <- subjectHits(ol)
    qhits <- queryHits(ol)
    local <- ranges(x)[qhits]
    bounds <- ranges(gr)[shits]

    ## location wrt start of individual list elements
    neg <- as.vector(strand(gr)[shits] == "-")
    local[!neg] <- shift(local[!neg], - start(bounds)[!neg])
    local[neg] <- IRanges(end(bounds)[neg] - end(local)[neg],
                          width = width(local)[neg])

    toInd <- togroup(to)[shits]
    ## location wrt start of list elements combined (e.g., transcript)
    cumsums <- .listCumsumShifted(width(to))
    txLevel <- shift(local, 1L + cumsums[shits])

    toInd <- togroup(to)[shits]
    if (elt.loc || elt.hits) {
        if (elt.loc && elt.hits) {
            mcols <- DataFrame(queryHits=qhits, subjectHits=toInd, 
                               eltLoc=local, eltHits=shits)
        } else {
        if (elt.loc)
            mcols <- DataFrame(queryHits=qhits, subjectHits=toInd,
                               eltLoc=shift(local, 1L))
        else
            mcols <- DataFrame(queryHits=qhits, subjectHits=toInd,
                               eltHits=shits)
        }
    } else {
        mcols <- DataFrame(queryHits=qhits, subjectHits=toInd)
    }
    GRanges(seqnames(gr)[shits], txLevel, strand = strand(gr[shits]), mcols)
  }
)
