### =========================================================================
### The seqinfo() (and releated) generic getters and setters
### -------------------------------------------------------------------------

setGeneric("seqinfo", function(x) standardGeneric("seqinfo"))

setGeneric("seqinfo<-", signature="x",
    function(x, old2new=NULL, value) standardGeneric("seqinfo<-")
)

setGeneric("seqnames", function(x) standardGeneric("seqnames"))

setGeneric("seqnames<-", signature="x",
    function(x, value) standardGeneric("seqnames<-")
)

setGeneric("seqlevels", function(x) standardGeneric("seqlevels"))

setMethod("seqlevels", "ANY", function(x) seqlevels(seqinfo(x)))

setGeneric("seqlevels<-", signature="x",
    function(x, value) standardGeneric("seqlevels<-")
)

setGeneric("seqlengths", function(x) standardGeneric("seqlengths"))

setMethod("seqlengths", "ANY", function(x) seqlengths(seqinfo(x)))

setGeneric("seqlengths<-", signature="x",
    function(x, value) standardGeneric("seqlengths<-")
)

setReplaceMethod("seqlengths", "ANY",
    function(x, value)
    {
        seqlengths(seqinfo(x)) <- value
        x
    }
)

setGeneric("isCircular", function(x) standardGeneric("isCircular"))

setMethod("isCircular", "ANY", function(x) isCircular(seqinfo(x)))

setGeneric("isCircular<-", signature="x",
    function(x, value) standardGeneric("isCircular<-")
)

setReplaceMethod("isCircular", "ANY",
    function(x, value)
    {
        isCircular(seqinfo(x)) <- value
        x
    }
)

setGeneric("isCircularWithKnownLength",
    function(x) standardGeneric("isCircularWithKnownLength")
)

setMethod("isCircularWithKnownLength", "ANY",
    function(x) isCircularWithKnownLength(seqinfo(x))
)

