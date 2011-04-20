### =========================================================================
### resolveHits methods
### -------------------------------------------------------------------------

setGeneric("resolveHits", signature = c("query", "subject"),
           function(query, subject,
                    type = c("any", "start", "end", "within", "equal"),
                    resolution = c("none", "divide", "uniqueDisjoint"),
                                   ignore.strand, ...)
           standardGeneric("resolveHits")
)

setMethod("resolveHits", c("GRangesList", "GRangesList"),
          function(query, subject,
              type = c("any", "start", "end", "within", "equal"),
              resolution = c("none", "divide", "uniqueDisjoint"), 
                             ignore.strand, ...)
{
    resolution <- match.arg(resolution)
    type <- match.arg(type)
    switch(resolution,
           ## FIXME: for 'none', what about 'clean' hits?
           none = rep.int(0L, length(unlist(subject))),
           divide = resDivide(query, subject, type, 
                              ignore.strand=ignore.strand),
           uniqueDisjoint = resUniqueDisjoint(query, subject, 
                              type, ignore.strand=ignore.strand))
})

resDivide <- function(query, subject, type, ignore.strand)
{
    fo <- findOverlaps(query, unlist(subject), type=type, 
                       ignore.strand=ignore.strand)
    qrle <- Rle(queryHits(fo))
    div <- rep.int((1/runLength(qrle)), runLength(qrle))
    idx <- unique(subjectHits(fo))
    counts <- lapply(split(div, subjectHits(fo)), sum)
    hits <- rep.int(0, length(unlist(subject)))
    hits[idx] <- unlist(counts) 
    hits
}

resUniqueDisjoint <- function(query, subject, type, ignore.strand)
{
    ## ud regions
    gr <- unlist(subject)
    d <- disjoin(gr)
    ud <- d[countOverlaps(d, gr, type=type, ignore.strand=ignore.strand) == 1]
    ## remove reads that hit multiple ud regions
    multihit <- countOverlaps(query, ud, type=type, ignore.strand=ignore.strand)
    fo <- findOverlaps(query[multihit == 1], ud, type=type,
        ignore.strand=ignore.strand)
    ## map ud regions back to original subjects
    backmap <- findOverlaps(ud[subjectHits(fo)], unlist(subject), type=type,
                       ignore.strand=ignore.strand)
    idx <- unique(subjectHits(backmap))
    counts <- lapply(split(subjectHits(backmap), subjectHits(backmap)), length)
    hits <- rep.int(0, length(unlist(subject)))
    hits[idx] <- unlist(counts)
    hits
}


