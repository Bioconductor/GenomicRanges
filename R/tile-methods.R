### =========================================================================
### "tile" methods
### -------------------------------------------------------------------------
###

setMethod("tile", "GenomicRanges", function(x, n, width) {
  sn <- seqnames(x)
  strand <- strand(x)
  x <- ranges(x)
  tiles <- callGeneric()
  gr <- GRanges(rep(sn, elementNROWS(tiles)), unlist(tiles),
                rep(strand, elementNROWS(tiles)))
  relist(gr, tiles)
})
