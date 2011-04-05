### =========================================================================
### resolveHits methods
### -------------------------------------------------------------------------

setGeneric("resolveHits", signature = c("query", "subject"),
           function(query, subject,
                    type = c("any", "start", "end", "within", "equal"),
                    resolution = c("none", "divide", "uniqueDisjoint",
                                   "coverage"), ...)
           standardGeneric("resolveHits")
)

setMethod("resolveHits", c("GRangesList", "GRangesList"),
          function(query, subject,
              type = c("any", "start", "end", "within", "equal"),
              resolution = c("none", "divide", "uniqueDisjoint", 
                            "coverage"), ...)
{
    resolution <- match.arg(resolution)
    type <- match.arg(type)
    switch(resolution,
           ## FIXME: for 'none', what about 'clean' hits?
           none = rep.int(0L, length(unlist(subject))),
           divide = resDivide(query, subject, type),
           uniqueDisjoint = resUniqueDisjoint(query, subject, type),
           coverage = resCoverage(query, subject, type)) 
})

resDivide <- function(query, subject, type)
{
    fo <- findOverlaps(query, unlist(subject), type=type)
    lst <- split(subjectHits(fo), queryHits(fo))
    len <- lapply(lst, length)
    div <- rep(1/unlist(len), unlist(len))
    hits <- rep.int(0, length(unlist(subject)))
    hits[unlist(lst)] <- div
    hits
}

resUniqueDisjoint <- function(query, subject, type)
{
    gr <- unlist(subject)
    d <- disjoin(gr)
    ud <- d[countOverlaps(d, gr) == 1]
    
    co <- countOverlaps(ud, query)
    fo <- findOverlaps(ud[co == 1], unlist(subject))
    el <- rep(elementLengths(subject), elementLengths(subject))
    h <- subjectHits(fo)[el[subjectHits(fo)] == 1]
    hits <- rep.int(0, length(unlist(subject)))
    hits[h] <- 1
    hits
}

.ud <- function(x)
{
    ## unique disjoint regions of GRangesList
    r <- GenomicRanges:::deconstructGRLintoGR(x)
    d <- disjoin(r)
    ud_p <- d[countOverlaps(d, r) == 1]
    rle <- seqnames(ud_p)
    s <- unlist(GenomicRanges:::reconstructGRLfromGR(ud_p, x))
    values(s) <- rle
    s
}

