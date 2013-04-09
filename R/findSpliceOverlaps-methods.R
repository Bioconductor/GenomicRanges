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

setMethod("findSpliceOverlaps", c("GAlignments", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
    callGeneric(grglist(query, order.as.in.query=TRUE), subject, 
                ignore.strand, ..., cds=cds)
})

setMethod("findSpliceOverlaps", c("GAlignmentPairs", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
### FIXME: order.as.in.query = FALSE needed for insertGaps(). If we
### really want to use insertGaps(), we need to make it robust to
### different orderings.

### FIXME:
### instead of relying on query.break column, maybe we should add a
### 'splice = .gaps(query)' argument to .findSpliceOverlaps that we
### set to introns(query) here. The downside is that a GRangesList
### derived from GAlignmentPairs will no longer work.
    callGeneric(grglist(query, order.as.in.query=FALSE), subject,
                ignore.strand, ..., cds=cds)
})

setMethod("findSpliceOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
    .findSpliceOverlaps(query, subject, ignore.strand, cds=cds)
})

.findSpliceOverlaps <- function(query, subject, ignore.strand=FALSE, cds=NULL)
{
    ## adjust strand based on 'XS'
    if (!is.null(xs <- mcols(query)$XS)) {
        strand <- ifelse(!is.na(xs), xs, "*")
        strand(query) <- relist(Rle(strand, elementLengths(query)), 
                                query)
    }
    ## NOTE: this misses reads completely within an intron, but this
    ## is intentional: a read is only assigned to a transcript if it
    ## hits an exon. Otherwise, it could be from another gene inside
    ## an intron (happens frequently).
    olap <- findOverlaps(query, subject, ignore.strand=ignore.strand)
    if (length(olap) == 0L)
        return(.result(olap))
    if (!is.null(cds)) {
        coding <- logical(length(olap))
        hits <- findOverlaps(query, cds, ignore.strand=ignore.strand)
        coding[queryHits(olap) %in% queryHits(hits)] <- TRUE
    } else {
        coding <- rep.int(NA, length(olap))
    }

    query <- query[queryHits(olap)]
    subject <- subject[subjectHits(olap)]
    splice <- .gaps(query)
    
    compatible <- .compatibleTranscription(query, subject, splice)
    unique <- .oneMatch(compatible, queryHits(olap))
    strandSpecific <- all(strand(query) != "*")
    mcols(olap) <- DataFrame(compatible, unique, coding, strandSpecific)
    olap 
}

## -------------------------------------------------------------------------
## Methods in Rsamtools :
## findSpliceOverlaps,character,ANY-method
## findSpliceOverlaps,BamFile,ANY-method
 
## -------------------------------------------------------------------------
## Methods in GenomicFeatures : 
## findSpliceOverlaps,GAlignments,TranscriptDb-method
## findSpliceOverlaps,GAlignmentPairs,TranscriptDb-method
## findSpliceOverlaps,GRangesList,TranscriptDb-method
