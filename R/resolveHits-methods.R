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
    if (type == "within") {
        fo <- findOverlaps(subject, query, type=type, ignore.strand=ignore.strand)
        hitFraction <- readValue * (1/table(queryHits(fo)))
        hitSplit <- split(hitFraction[queryHits(fo)], subjectHits(fo))
    } else { 
        fo <- findOverlaps(query, subject, type=type, ignore.strand=ignore.strand)
        hitFraction <- readValue * (1/table(subjectHits(fo)))
        hitSplit <- split(hitFraction[subjectHits(fo)], queryHits(fo))
    }
    hits <- double(length(query))
    hits[as.numeric(names(hitSplit))] <- unlist(lapply(hitSplit, sum))
    hits
}

resUniqueDisjoint <- function(query, subject, readValue, type, ignore.strand)
{
    ## ud regions
    d <- disjoin(query)
    #ud <- d[countOverlaps(d, query, type=type, ignore.strand=ignore.strand) == 1]
    ud <- d[countOverlaps(d, query, type="any", ignore.strand=ignore.strand) == 1]
    ## ignore reads that hit multiple ud regions
    multihit <- countOverlaps(subject, ud, type="any", ignore.strand=ignore.strand)

    if (type == "within") {
        fo <- findOverlaps(subject[multihit == 1], ud, type=type,
            ignore.strand=ignore.strand)
        ## map ud regions back to original query 
        backmap <- findOverlaps(ud[subjectHits(fo)], query, type="any",
            ignore.strand=ignore.strand)
        hitSplit <- split(readValue[multihit == 1][queryHits(fo)], 
            subjectHits(backmap))
    } else {
        fo <- findOverlaps(ud, subject[multihit == 1], type=type,
            ignore.strand=ignore.strand)
        ## map ud regions back to original query 
        backmap <- findOverlaps(ud[queryHits(fo)], query, type="any",
            ignore.strand=ignore.strand)
        hitSplit <- split(readValue[multihit == 1][subjectHits(fo)], 
            subjectHits(backmap))
    }
    hits <- double(length(query))
    hits[as.numeric(names(hitSplit))] <- unlist(lapply(hitSplit, sum))
    hits
}

