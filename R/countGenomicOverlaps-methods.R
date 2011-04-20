### =========================================================================
### countGenomicOverlaps methods
### -------------------------------------------------------------------------

setGeneric("countGenomicOverlaps", signature = c("query", "subject"), 
            function(query, subject, 
                     type = c("any", "start", "end", "within", "equal"),
                     resolution = c("none", "divide", "uniqueDisjoint"),
                     junctionwt = 0.5,
                     ignore.strand = FALSE, ...)
            standardGeneric("countGenomicOverlaps")
)

setMethod("countGenomicOverlaps", c("GRangesList", "GenomicRanges"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             junctionwt = 0.5,
             ignore.strand = FALSE, ...)
{
    listSubject <- split(subject, seq_len(length(subject))) 
    callGeneric(query, listSubject, 
                type = type, resolution = resolution, 
                junctionwt = junctionwt, ignore.strand = ignore.strand, ...)
})

setMethod("countGenomicOverlaps", c("GenomicRanges", "GRangesList"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             junctionwt = 0.5,
             ignore.strand = FALSE, ...)
{
    listQuery <- split(query, seq_len(length(query))) 
    callGeneric(listQuery, subject, 
                type = type, resolution = resolution, 
                junctionwt = junctionwt, ignore.strand = ignore.strand, ...)
})

setMethod("countGenomicOverlaps", c("GenomicRanges", "GenomicRanges"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             junctionwt = 0.5,
             ignore.strand = FALSE, ...)
{
    listSubject <- split(subject, seq_len(length(subject))) 
    listQuery <- split(query, seq_len(length(query))) 
    callGeneric(listQuery, listSubject, 
                type = type, resolution = resolution, 
                junctionwt = junctionwt, ignore.strand = ignore.strand, ...)
})

setMethod("countGenomicOverlaps", c("GappedAlignments", "GRangesList"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             junctionwt = 0.5,
             ignore.strand = FALSE, ...)
{
    listQuery <- as(query, "GRangesList") 
    callGeneric(listQuery, subject, 
                type = type, resolution = resolution, 
                junctionwt = junctionwt, ignore.strand = ignore.strand, ...)
})

setMethod("countGenomicOverlaps", c("GappedAlignments", "GenomicRanges"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             junctionwt = 0.5,
             ignore.strand = FALSE, ...)
{
    listSubject <- split(subject, seq_len(length(subject))) 
    listQuery <- as(query, "GRangesList") 
    callGeneric(listQuery, listSubject, 
                type = type, resolution = resolution, 
                junctionwt = junctionwt, ignore.strand = ignore.strand, ...)
})

setMethod("countGenomicOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             junctionwt = 0.5,
             ignore.strand = FALSE, ...)
{
    resolution <- match.arg(resolution)
    type <- match.arg(type)
    if (ignore.strand)
        strand(subject@unlistData) <- "*"
    if (type == "within" && resolution == "uniqueDisjoint")
        stop("resolution `uniqueDisjoint' with type `within'",
             "is not logical")
    co <- countOverlaps(query, unlist(subject), type=type, 
                        ignore.strand=ignore.strand)
    if (!any(co)) {
        values(subject@unlistData)[["hits"]] <- 
          rep.int(0L, length(unlist(subject))) 
        return(subject)
    }

    ## split reads
    if (resolution == "none") {
        split <- rep.int(0, length(unlist(subject))) 
        query <- query[elementLengths(query) <= 1]
    } else {
        if (any(elementLengths(query) > 1)) { 
            sr <- query[elementLengths(query) > 1]
            co_split <- countOverlaps(unlist(sr), unlist(subject), type=type,
                                ignore.strand=ignore.strand)
            if (any(co_split > 2)) {
                warning("split reads that hit > 2 exons have been ",
                        "detected and will be dropped")
                query <- query[co_split < 3] 
                sr <- query[elementLengths(query) > 1]
            } 
            ## exon hit by both read seqments gets wt = 1 
            ## exons hit by one read seqments each get wt = junctionwt 
            fo <- findOverlaps(unlist(sr), unlist(subject), type=type,
                               ignore.strand=ignore.strand)
            subRle <- Rle(subjectHits(fo))
            wt <- runLength(subRle)
            wt <- ifelse(wt == 1, junctionwt, 1)
            split <- rep.int(0, length(unlist(subject)))
            split[runValue(subRle)] <- wt 
            query <- query[elementLengths(query) <= 1]
        } else split <- rep.int(0, length(unlist(subject)))
    } 
 
    co <- countOverlaps(query, unlist(subject), type=type, 
                        ignore.strand=ignore.strand)
    ## read hit one subject
    if (any(co == 1)) {
        clean <- countOverlaps(unlist(subject), query[co == 1],
                               ignore.strand=ignore.strand)
    } else {
        clean <- rep.int(0, length(unlist(subject)))
    }
    ## read hit multiple subjects
    if (any(co > 1)) {
        resolved <- resolveHits(query[co > 1], subject, type=type, 
                                resolution=resolution, 
                                ignore.strand=ignore.strand) 
    } else {
        resolved <- rep.int(0, length(unlist(subject)))
    }
    values(subject@unlistData)[["hits"]] <- split + clean + resolved 
    subject
})
