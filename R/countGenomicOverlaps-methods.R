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
    listSubject <- split(subject, seq_len(length(subject))) 
    listQuery <- split(query, seq_len(length(query))) 
    callGeneric(listQuery, listSubject, type = type, resolution = resolution, 
        ignore.strand = ignore.strand, ...)
})

setMethod("countGenomicOverlaps", c("GappedAlignments", "GRangesList"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, splitreads = TRUE, ...)
{
    listQuery <- as(query, "GRangesList") 
    callGeneric(listQuery, subject, type = type, resolution = resolution, 
        ignore.strand = ignore.strand, splitreads = TRUE, ...)
})

setMethod("countGenomicOverlaps", c("GappedAlignments", "GenomicRanges"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, splitreads = TRUE, ...)
{
    listSubject <- split(subject, seq_len(length(subject))) 
    listQuery <- as(query, "GRangesList") 
    callGeneric(listQuery, listSubject, type = type, resolution = resolution, 
        ignore.strand = ignore.strand, splitreads = TRUE, ...)
})

setMethod("countGenomicOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, 
             type = c("any", "start", "end", "within", "equal"),
             resolution = c("none", "divide", "uniqueDisjoint"), 
             ignore.strand = FALSE, splitreads = TRUE, ...)
{
    resolution <- match.arg(resolution)
    type <- match.arg(type)
    counts <- .countGenomicOverlaps(query, subject, 
        type = type, resolution = resolution, 
        ignore.strand = ignore.strand, splitreads = splitreads, ...)
    if (length(metadata(query)) > 0)
        colData <- DataFrame(metaData = metadata(query))
    else
        colData <- DataFrame(metaData = character(1))
    SummarizedExperiment(assays=SimpleList(counts=counts),
                         ## FIXME : remove unlist when rowData can take GRList 
                         rowData=unlist(subject), colData=colData)
})

.countGenomicOverlaps <- function(query, subject, type, resolution,
    ignore.strand, splitreads)
{
    if (ignore.strand)
        strand(subject@unlistData) <- "*"
    if (type == "within" && resolution == "uniqueDisjoint")
        stop("resolution `uniqueDisjoint' with type `within'",
             "is not logical")

    ## value each read or read fragment contributes
    if (splitreads == TRUE) {
        readValue <- rep.int(1/elementLengths(query), elementLengths(query))
        uquery <- unlist(query, use.names=FALSE)
    } else {
        query <- query[elementLengths(query) == 1]
        readValue <- rep.int(1, length(query))
        uquery <- unlist(query, use.names=FALSE)
    }

    usubject <- unlist(subject, use.names=FALSE)
    co <- countOverlaps(uquery, usubject, 
        type=type, ignore.strand=ignore.strand)
    if (!any(co)) {
        values(subject@unlistData)[["hits"]] <- 
          integer(length(usubject)) 
        warning("no overlaps detected")
        return(subject)
    }

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

    matrix(clean + resolved, ncol=1) 
}
