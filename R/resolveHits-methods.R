### =========================================================================
### resolveHits methods
### -------------------------------------------------------------------------

setGeneric("resolveHits", signature = c("query", "subject"),  # not exported
           function(query, subject, readValue,
                    type = c("any", "start", "end", "within", "equal"),
                    resolution = c("divide", "uniqueDisjoint"),
                                   ignore.strand, ...)
           standardGeneric("resolveHits")
)

setMethod("resolveHits", c("GenomicRanges", "GenomicRanges"),  # not exported
          function(query, subject, readValue,
              type = c("any", "start", "end", "within", "equal"),
              resolution = c("divide", "uniqueDisjoint"), 
                             ignore.strand, ...)
{
    .Deprecated("summarizeOverlaps")
})
