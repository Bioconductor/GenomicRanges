### =========================================================================
### RangedData/RangesList implementation of the GenomicRanges API
### -------------------------------------------------------------------------
###

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqinfo
###

setMethod("seqinfo", "List", function(x) {
  si <- metadata(x)$seqinfo
  if (is.null(si)) {
    sn <- names(x)
    if (is.null(sn))
      sn <- as.character(seq(length(x)))
    si <- Seqinfo(sn)
  }
  si  
})

setMethod("seqinfo", "RangesList", function(x) {
  si <- callNextMethod()
  if (!is.null(universe(x)))
    genome(si) <- universe(x)
  si
})

### FIXME: needs sanity checks
setReplaceMethod("seqinfo", "List",
                 function(x, value) {
                   metadata(x)$seqinfo <- value
                   x
                 })

setMethod("seqinfo", "RangedData", function(x) seqinfo(ranges(x)))
setReplaceMethod("seqinfo", "RangedData",
                 function(x, value) {
                   seqinfo(ranges(x)) <- value
                   x
                 })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqnames
###

setMethod("seqnames", "RangesList", function(x) {
  if (is.null(names(x)))
    NULL
  else seqsplit(Rle(space(x)), rep(names(x), elementLengths(x)))
})

setMethod("seqnames", "RangedData", function(x) space(x))

