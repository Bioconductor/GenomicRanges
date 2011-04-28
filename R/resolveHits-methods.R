### =========================================================================
### resolveHits methods
### -------------------------------------------------------------------------

setGeneric("resolveHits", signature = c("query", "subject"),
           function(query, subject, readValue,
                    type = c("any", "start", "end", "within", "equal"),
                    resolution = c("divide", "uniqueDisjoint"),
                                   ignore.strand, ...)
           standardGeneric("resolveHits")
)

setMethod("resolveHits", c("GenomicRanges", "GenomicRanges"),
          function(query, subject, readValue,
              type = c("any", "start", "end", "within", "equal"),
              resolution = c("divide", "uniqueDisjoint"), 
                             ignore.strand, ...)
{
    resolution <- match.arg(resolution)
    type <- match.arg(type)
    switch(resolution,
           divide = resDivide(query, subject, readValue, type, 
               ignore.strand=ignore.strand),
           uniqueDisjoint = resUniqueDisjoint(query, subject, readValue, 
               type, ignore.strand=ignore.strand))
})

resDivide <- function(query, subject, readValue, type, ignore.strand)
{
    fo <- findOverlaps(query, subject, type=type, ignore.strand=ignore.strand)
    hitFraction <- readValue * (1/table(queryHits(fo)))
    hitSplit <- split(hitFraction[queryHits(fo)], subjectHits(fo))
    hits <- double(length(subject))
    hits[as.numeric(names(hitSplit))] <- unlist(lapply(hitSplit, sum))
    hits
}

resUniqueDisjoint <- function(query, subject, readValue, type, ignore.strand)
{
    ## ud regions
    d <- disjoin(subject)
    ud <- d[countOverlaps(d, subject, type=type, ignore.strand=ignore.strand) == 1]
    ## ignore reads that hit multiple ud regions
    multihit <- countOverlaps(query, ud, type=type, ignore.strand=ignore.strand)
    fo <- findOverlaps(query[multihit == 1], ud, type=type,
        ignore.strand=ignore.strand)
    ## map ud regions back to original subjects
    backmap <- findOverlaps(ud[subjectHits(fo)], subject, type=type,
        ignore.strand=ignore.strand)

    hitSplit <- split(readValue[multihit == 1][queryHits(fo)], 
        subjectHits(backmap))
    hits <- double(length(subject))
    hits[as.numeric(names(hitSplit))] <- unlist(lapply(hitSplit, sum))
    hits
}

