### =========================================================================
### 'mapToGenome' and 'mapToTranscript' methods
### -------------------------------------------------------------------------
###

### Generics.

setGeneric("mapToGenome", signature=c("x", "alignment"),
    function(x, alignment, ...) standardGeneric("mapToGenome")
)

setGeneric("pmapToGenome", signature=c("x", "alignment"),
    function(x, alignment, ...) standardGeneric("pmapToGenome")
)

setGeneric("mapToTranscript", signature=c("x", "alignment"),
    function(x, alignment, ...) standardGeneric("mapToTranscript")
)

setGeneric("pmapToTranscript", signature=c("x", "alignment"),
    function(x, alignment, ...) standardGeneric("pmapToTranscript")
)

