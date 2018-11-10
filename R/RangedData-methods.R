### =========================================================================
### RangedData/IntegerRangesList implementation of the GenomicRanges API
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
      sn <- as.character(seq_len(length(x)))
    si <- Seqinfo(unique(sn))
  }
  si  
})

### FIXME: needs sanity checks
setReplaceMethod("seqinfo", "List",
                 function(x, value) {
                   metadata(x)$seqinfo <- value
                   x
                 })

setMethod("seqinfo", "RangedData",
  function(x)
  {
    .Deprecated(msg=wmsg(IRanges:::RangedData_is_deprecated_msg))
    seqinfo(ranges(x))
  }
)

setReplaceMethod("seqinfo", "RangedData",
                 function(x, value) {
                   what <- "seqinfo() setter for RangedData objects"
                   .Defunct(msg=wmsg(IRanges:::RangedData_method_is_defunct_msg(what)))
                   seqinfo(ranges(x)) <- value
                   x
                 })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqnames
###

setMethod("seqnames", "RangedData",
  function(x)
  {
    what <- "seqnames() getter for RangedData objects"
    .Defunct(msg=wmsg(IRanges:::RangedData_method_is_defunct_msg(what)))
    space(x)
  }
)

