### The seqnames() generic and methods.

setGeneric("seqnames", function(x) standardGeneric("seqnames"))
setGeneric("seqnames<-", function(x, value) standardGeneric("seqnames<-"))

### The seqlengths() generic and methods.

setGeneric("seqlengths", function(x) standardGeneric("seqlengths"))
setGeneric("seqlengths<-", function(x, value) standardGeneric("seqlengths<-"))
