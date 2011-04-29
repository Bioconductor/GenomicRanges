### =========================================================================
### countGenomicOverlaps methods
### -------------------------------------------------------------------------

setGeneric("countGenomicOverlaps", signature = c("query", "subject"), 
            function(query, subject, 
                     type = c("any", "start", "end", "within", "equal"),
                     resolution = c("none", "divide", "uniqueDisjoint"),
                     ignore.strand = FALSE, ...)
            standardGeneric("countGenomicOverlaps")
)

setMethod("countGenomicOverlaps", c("GRangesList", "GenomicRanges"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, ...)
{
    listSubject <- split(subject, seq_len(length(subject))) 
    ans <- callGeneric(query, listSubject, 
        type = type, resolution = resolution, 
        ignore.strand = ignore.strand, ...)
    unlist(ans)
})

setMethod("countGenomicOverlaps", c("GenomicRanges", "GRangesList"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, ...)
{
    listQuery <- split(query, seq_len(length(query))) 
    callGeneric(listQuery, subject, 
                type = type, resolution = resolution, 
                ignore.strand = ignore.strand, ...)
})

setMethod("countGenomicOverlaps", c("GenomicRanges", "GenomicRanges"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, ...)
{
    listSubject <- split(subject, seq_len(length(subject))) 
    listQuery <- split(query, seq_len(length(query))) 
    ans <- callGeneric(listQuery, listSubject, 
        type = type, resolution = resolution, 
        ignore.strand = ignore.strand, ...)
    unlist(ans)
})

setMethod("countGenomicOverlaps", c("GappedAlignments", "GRangesList"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, ...)
{
    listQuery <- as(query, "GRangesList") 
    callGeneric(listQuery, subject, 
                type = type, resolution = resolution, 
                ignore.strand = ignore.strand, ...)
})

setMethod("countGenomicOverlaps", c("GappedAlignments", "GenomicRanges"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, ...)
{
    listSubject <- split(subject, seq_len(length(subject))) 
    listQuery <- as(query, "GRangesList") 
    ans <- callGeneric(listQuery, listSubject, 
        type = type, resolution = resolution, 
        ignore.strand = ignore.strand, ...)
    unlist(ans)
})

setMethod("countGenomicOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, ...)
{
    resolution <- match.arg(resolution)
    type <- match.arg(type)
    if (ignore.strand)
        strand(subject@unlistData) <- "*"
    if (type == "within" && resolution == "uniqueDisjoint")
        stop("resolution `uniqueDisjoint' with type `within'",
             "is not logical")

    uquery <- unlist(query, use.names=FALSE)
    usubject <- unlist(subject, use.names=FALSE)

    co <- countOverlaps(uquery, usubject, 
        type=type, ignore.strand=ignore.strand)
    if (!any(co)) {
        values(subject@unlistData)[["hits"]] <- 
          integer(length(usubject)) 
        return(subject)
    }

    ## define the hit value each read contributes
    readValue <- rep.int(1/elementLengths(query), elementLengths(query))

    ## read hit one subject
    if (any(co == 1)) {
    fo <- findOverlaps(uquery[co == 1], usubject, type=type,
        ignore.strand=ignore.strand)
    cleanSplit <- split(readValue[co == 1][queryHits(fo)], subjectHits(fo))
    clean <- double(length(usubject))
    clean[as.numeric(names(cleanSplit))] <- unlist(lapply(cleanSplit, sum)) 
    } else {
        clean <- double(length(usubject))
    }

    ## read hit multiple subjects
    if (any(co > 1) && resolution != "none") {
        resolved <- resolveHits(uquery[co > 1], usubject, readValue[co > 1],
            type=type, resolution=resolution, ignore.strand=ignore.strand) 
    } else {
        resolved <- double(length(usubject))
    }
    values(subject@unlistData)[["hits"]] <- clean + resolved 
    subject
})
