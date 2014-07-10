### =========================================================================
### "map" methods
### -------------------------------------------------------------------------
###


.listCumsum <- function(x) {
  x_unlisted <- unlist(x, use.names=FALSE)
  x_cumsum <- cumsum(x_unlisted)
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
  neg <- strand(gr)[start(part)] == "-"
  ord <- S4Vectors:::mseq(ifelse(neg, end(part), start(part)),
                          ifelse(neg, start(part), end(part)))
  relist(gr[ord], x)
}

setMethod("map", c("GenomicRanges", "GRangesList"), function(from, to) {
  ## make sure 'to' is properly sorted by strand
  to <- .orderElementsByTranscription(to)

  ## find overlaps
  gr <- unlist(to, use.names=FALSE)
  ol <- findOverlaps(from, gr, type = "within")
  shits <- subjectHits(ol)
  qhits <- queryHits(ol)
  local <- ranges(from)[qhits]
  bounds <- ranges(gr)[shits]

  ## location wrt start of coding region 
  neg <- as.vector(strand(gr)[shits] == "-")
  local[!neg] <- shift(local[!neg], - start(bounds)[!neg])
  local[neg] <- IRanges(end(bounds)[neg] - end(local)[neg],
                        width = width(local)[neg])

  ## location wrt transcript
  cumsums <- .listCumsumShifted(width(to))
  local <- shift(local, 1L + cumsums[shits])

  toInd <- togroup(to)[shits]
  matching <- new("Hits",
                  queryHits = qhits, subjectHits = toInd,
                  queryLength = length(from), subjectLength = length(to))

  new("GRangesMapping", hits = matching,
      granges = GRanges(seqnames(gr)[shits], local, 
                        strand = strand(gr[shits])))
})

