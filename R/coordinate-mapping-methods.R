### =========================================================================
### 'mapToGenome' and 'mapToTranscript' methods
### -------------------------------------------------------------------------
###

### Generics.

setGeneric("mapToGenome", signature=c("from", "to"),
    function(from, to, ...) standardGeneric("mapToGenome")
)

setGeneric("pmapToGenome", signature=c("from", "to"),
    function(from, to, ...) standardGeneric("pmapToGenome")
)

setGeneric("mapToTranscript", signature=c("from", "to"),
    function(from, to, ...) standardGeneric("mapToTranscript")
)

setGeneric("pmapToTranscript", signature=c("from", "to"),
    function(from, to, ...) standardGeneric("pmapToTranscript")
)

