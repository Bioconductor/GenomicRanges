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
    hits <- rep.int(0, length(unlist(subject)))
    hits[subjectHits(fo)] <- div
    hits
}

resUniqueDisjoint <- function(query, subject, type, ignore.strand)
{
    gr <- unlist(subject)
    d <- disjoin(gr)
    ud <- d[countOverlaps(d, gr, type=type, ignore.strand=ignore.strand) == 1]
 
    co <- countOverlaps(ud, query, type=type, ignore.strand=ignore.strand)
    fo <- findOverlaps(ud[co == 1], unlist(subject), type=type,
                       ignore.strand=ignore.strand)
    el <- rep(elementLengths(subject), elementLengths(subject))
    h <- subjectHits(fo)[el[subjectHits(fo)] == 1]
    hits <- rep.int(0, length(unlist(subject)))
    hits[h] <- 1
    hits
}


