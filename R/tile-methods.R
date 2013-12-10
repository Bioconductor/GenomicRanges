### =========================================================================
### "tile" methods
### -------------------------------------------------------------------------
###

setMethod("tile", "GenomicRanges", function(x, n, width) {
  sn <- seqnames(x)
  strand <- strand(x)
  x <- ranges(x)
  tiles <- callGeneric()
  gr <- GRanges(rep(sn, elementLengths(tiles)), unlist(tiles),
                rep(strand, elementLengths(tiles)))
  relist(gr, tiles)
})
