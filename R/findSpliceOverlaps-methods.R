### =========================================================================
### findSpliceOverlaps methods 
### -------------------------------------------------------------------------

setGeneric("findSpliceOverlaps", signature=c("query", "subject"),
    function(query, subject, ignore.strand=FALSE, ...)
{
    standardGeneric("findSpliceOverlaps")
})

## -------------------------------------------------------------------------
## subject as GRangesList 
##

setMethod("findSpliceOverlaps", c("GappedAlignments", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
    callGeneric(grglist(query, order.as.in.query=TRUE), subject, 
                ignore.strand, ..., cds=cds)
})

setMethod("findSpliceOverlaps", c("GappedAlignmentPairs", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
    .findSpliceOverlaps(query, subject, ignore.strand, pairedEnd=TRUE, cds=cds)
})

setMethod("findSpliceOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
    .findSpliceOverlaps(query, subject, ignore.strand, cds=cds)
})

.findSpliceOverlaps <- function(query, subject, ignore.strand=FALSE,
                                pairedEnd=FALSE, cds=NULL)
{
    ## FIXME : this overlap misses a read competely within an intron
    olap <- findOverlaps(query, subject, ignore.strand=ignore.strand)
    if (length(olap) == 0L)
        return(.result(olap))
    if (!is.null(cds)) {
        coding <- logical(length(olap))
        hits <- findOverlaps(subject, cds, ignore.strand=ignore.strand)
        coding[subjectHits(olap) %in% queryHits(hits)] <- TRUE
    } else {
        coding <- rep.int(NA, length(olap))
    }
    intron <- .gaps(subject)
    query <- query[queryHits(olap)]
    subject <- subject[subjectHits(olap)]
    intron <- intron[subjectHits(olap)]
    if (pairedEnd) {
        query <- grglist(query, order.as.in.query=TRUE)
    } else { 
        splice <- .gaps(query)
    }
    compatible <- .compatibleTranscription(query, subject, splice)
    unique <- .oneMatch(compatible, queryHits(olap))
    strandSpecific <- all(strand(query) != "*")[queryHits(olap)]
    values(olap) <- DataFrame(compatible, unique, coding, strandSpecific)
    olap 
}

## -------------------------------------------------------------------------
## Methods in Rsamtools :
## findSpliceOverlaps,character,ANY-method
## findSpliceOverlaps,BamFile,ANY-method
 
## -------------------------------------------------------------------------
## Methods in GenomicFeatures : 
## findSpliceOverlaps,GappedAlignments,TranscriptDb-method
## findSpliceOverlaps,GappedAlignmentPairs,TranscriptDb-method
## findSpliceOverlaps,GRangesList,TranscriptDb-method
