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

setMethod("countGenomicOverlaps", c("GRangesList", "GRanges"),
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

setMethod("countGenomicOverlaps", c("GRanges", "GRangesList"),
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

setMethod("countGenomicOverlaps", c("GRanges", "GRanges"),
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

setMethod("countGenomicOverlaps", c("GappedAlignments", "GRanges"),
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
    if (!ignore.strand)
        strand(subject@unlistData) <- "*"
    if (type == "within" && resolution == "uniqueDisjoint")
        stop("It is not logical to use resolution `uniqueDisjoint` with ",
             "type `within`.")
    co <- countOverlaps(query, unlist(subject), type=type)
    if (!any(co)) {
        values(subject@unlistData)[["hits"]] <- 
          rep.int(0L, length(unlist(subject))) 
        return(subject)
    }

    ## split reads
    if (resolution == "none") {
        split <- rep.int(0, length(unlist(subject))) 
        query <- query[!(elementLengths(query) > 1)]
    } else {
        if (any(elementLengths(query) > 1)) { 
            sr <- query[elementLengths(query) > 1]
            co <- countOverlaps(unlist(sr), unlist(subject), type=type)
            if (any(co > 2)) {
                warning("Split reads that hit > 2 exons have been ",
                        "detected and will be dropped.")
                query <- query[co < 3] 
                sr <- query[elementLengths(query) > 1]
            } 
            ## exon hit by both read seqments gets wt = 1 
            ## exons hit by one read seqments each get wt = junctionwt 
            fo <- findOverlaps(unlist(sr), unlist(subject), type=type)
            subRle <- Rle(subjectHits(fo))
            wt <- runLength(subRle)
            wt <- ifelse(wt == 1, junctionwt, 1)
            split <- rep.int(0, length(unlist(subject)))
            split[runValue(subRle)] <- wt 
            query <- query[!(elementLengths(query) > 1)]
        } else split <- rep.int(0, length(unlist(subject)))
    } 
 
    co <- countOverlaps(query, unlist(subject), type=type)
    ## read hit one subject
    if (any(co == 1)) {
        clean <- countOverlaps(unlist(subject), query[co == 1])
    } else {
        clean <- rep.int(0, length(unlist(subject)))
    }
    ## read hit multiple subjects
    if (any(co > 1)) {
        resolved <- resolveHits(query[co > 1], subject, type=type, 
                                resolution=resolution) 
    } else {
        resolved <- rep.int(0, length(unlist(subject)))
    }
    values(subject@unlistData)[["hits"]] <- split + clean + resolved 
    subject
})
