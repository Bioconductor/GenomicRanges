### =========================================================================
### "tile" and "slidingWindows" methods
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

setMethod("slidingWindows", "GenomicRanges", function(x, width, step = 1L) {
    windows <- slidingWindows(ranges(x), width, step)
    gr <- rep(granges(x), lengths(windows))
    ranges(gr) <- unlist(windows)
    relist(gr, windows)
})
