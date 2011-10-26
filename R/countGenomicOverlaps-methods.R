### =========================================================================
### countGenomicOverlaps methods
### -------------------------------------------------------------------------

setGeneric("countGenomicOverlaps", signature = c("query", "subject"), 
            function(query, subject, 
                     type = c("any", "start", "end", "within", "equal"),
                     resolution = c("none", "divide", "uniqueDisjoint"),
                     ignore.strand = FALSE, splitreads = TRUE, ...)
            standardGeneric("countGenomicOverlaps")
)

setMethod("countGenomicOverlaps", c("GRangesList", "GenomicRanges"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, splitreads = TRUE, ...)
{
    msg <- c(" countGenomicOverlaps is deprecated.\n Please use ",
             "summarizeOverlaps instead.")
    .Deprecated(msg=msg)
    listSubject <- split(subject, seq_len(length(subject))) 
    callGeneric(query, listSubject, type = type, resolution = resolution, 
        ignore.strand = ignore.strand, splitreads = TRUE, ...)
})

setMethod("countGenomicOverlaps", c("GenomicRanges", "GRangesList"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, splitreads = TRUE, ...)
{
    msg <- c(" countGenomicOverlaps is deprecated.\n Please use ",
             "summarizeOverlaps instead.")
    .Deprecated(msg=msg)
    listQuery <- split(query, seq_len(length(query))) 
    callGeneric(listQuery, subject, type = type, resolution = resolution, 
        ignore.strand = ignore.strand, ...)
})

setMethod("countGenomicOverlaps", c("GenomicRanges", "GenomicRanges"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, splitreads = TRUE, ...)
{
    msg <- c(" countGenomicOverlaps is deprecated.\n Please use ",
             "summarizeOverlaps instead.")
    .Deprecated(msg=msg)
    listSubject <- split(subject, seq_len(length(subject))) 
    listQuery <- split(query, seq_len(length(query))) 
    callGeneric(listQuery, listSubject, type = type, resolution = resolution, 
        ignore.strand = ignore.strand, ...)
})

setMethod("countGenomicOverlaps", c("GRangesList", "GappedAlignments"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, splitreads = TRUE, ...)
{
    msg <- c(" countGenomicOverlaps is deprecated.\n Please use ",
             "summarizeOverlaps instead.")
    .Deprecated(msg=msg)
    listSubject <- as(subject, "GRangesList") 
    callGeneric(query, listSubject, type = type, resolution = resolution, 
        ignore.strand = ignore.strand, splitreads = TRUE, ...)
})

setMethod("countGenomicOverlaps", c("GenomicRanges", "GappedAlignments"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, splitreads = TRUE, ...)
{
    msg <- c(" countGenomicOverlaps is deprecated.\n Please use ",
             "summarizeOverlaps instead.")
    .Deprecated(msg=msg)
    listQuery <- split(query, seq_len(length(query))) 
    listSubject <- as(subject, "GRangesList") 
    callGeneric(listQuery, listSubject, type = type, resolution = resolution, 
        ignore.strand = ignore.strand, splitreads = TRUE, ...)
})

setMethod("countGenomicOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, splitreads = TRUE, ...)
{
    msg <- c(" countGenomicOverlaps is deprecated.\n Please use ",
             "summarizeOverlaps instead.")
    .Deprecated(msg=msg)
    resolution <- match.arg(resolution)
    type <- match.arg(type)
    counts <- .countGenomicOverlaps(query, subject, 
        type = type, resolution = resolution, 
        ignore.strand = ignore.strand, splitreads = splitreads, ...)
})

.countGenomicOverlaps <- function(query, subject, type, resolution,
    ignore.strand, splitreads)
{
    if (ignore.strand)
        strand(query@unlistData) <- "*"
    if (type == "within" && resolution == "uniqueDisjoint")
        stop("resolution `uniqueDisjoint' with type `within'",
             "is not logical")

    ## value each read or read fragment contributes
    if (splitreads == TRUE) {
        readValue <- rep.int(1/elementLengths(subject), elementLengths(subject))
        usubject <- unlist(subject, use.names=FALSE)
    } else {
        subject <- subject[elementLengths(subject) == 1]
        readValue <- rep.int(1, length(subject))
        usubject <- unlist(subject, use.names=FALSE)
    }

    uquery <- unlist(query, use.names=FALSE)
    co <- countOverlaps(usubject, uquery, 
        type=type, ignore.strand=ignore.strand)
    if (!any(co)) {
        warning("no overlaps detected")
        return(integer(length(uquery)))
    }

    ## read hit one subject
    if (any(co == 1)) {
        ## type="within" handle separately
        ## findOverlaps type="within" is directional (i.e., query within subject)
        if (type == "within") {
            fo <- findOverlaps(usubject[co == 1], uquery, type=type,
                ignore.strand=ignore.strand)
            cleanSplit <- split(readValue[co == 1][queryHits(fo)], subjectHits(fo))
        } else {
            fo <- findOverlaps(uquery, usubject[co == 1], type=type,
                ignore.strand=ignore.strand)
            cleanSplit <- split(readValue[co == 1][subjectHits(fo)], queryHits(fo))
        }
        clean <- double(length(uquery))
        clean[as.numeric(names(cleanSplit))] <- unlist(lapply(cleanSplit, sum)) 
    } else {
        clean <- double(length(uquery))
    }

    ## read hit multiple subjects
    if (any(co > 1) && resolution != "none") {
        resolved <- resolveHits(uquery, usubject[co > 1], readValue[co > 1],
            type=type, resolution=resolution, ignore.strand=ignore.strand) 
    } else {
        resolved <- double(length(uquery))
    }

    matrix(clean + resolved, ncol=1) 
}
