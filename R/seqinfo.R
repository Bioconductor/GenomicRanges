### =========================================================================
### The seqinfo() (and releated) generic getters and setters
### -------------------------------------------------------------------------

setGeneric("seqinfo", function(x) standardGeneric("seqinfo"))

setGeneric("seqinfo<-", function(x, value) standardGeneric("seqinfo<-"))

setGeneric("seqnames", function(x) standardGeneric("seqnames"))

setGeneric("seqnames<-", function(x, value) standardGeneric("seqnames<-"))

setGeneric("seqlengths", function(x) standardGeneric("seqlengths"))

setMethod("seqlengths", "ANY", function(x) seqlengths(seqinfo(x)))

setGeneric("seqlengths<-", function(x, value) standardGeneric("seqlengths<-"))

setGeneric("isCircular", function(x) standardGeneric("isCircular"))

setMethod("isCircular", "ANY", function(x) isCircular(seqinfo(x)))

setGeneric("isCircular<-", function(x, value) standardGeneric("isCircular<-"))

setGeneric("isCircularWithKnownLength",
    function(x) standardGeneric("isCircularWithKnownLength")
)

setMethod("isCircularWithKnownLength", "ANY",
    function(x) ((isCircular(x) %in% TRUE) & !is.na(seqlengths(x)))
)

