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
    callGeneric(grglist(query, order.as.in.query=TRUE), subject,
                ignore.strand, ..., cds=cds)
})

setMethod("findSpliceOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
    .findSpliceOverlaps(query, subject, ignore.strand, cds=cds)
})

.findSpliceOverlaps <- function(query, subject, ignore.strand=FALSE, cds=NULL)
{
    ## FIXME : ignore.strand?
    olap <- findOverlaps(query, subject, ignore.strand=ignore.strand)
    if (length(olap) == 0L)
        return(.result(olap))
    query <- query[queryHits(olap)]
    subject <- subject[subjectHits(olap)]

    if (!is.null(cds)) {
        coding <- rep.int(FALSE, length(olap))
        hits <- findOverlaps(subject, cds, ignore.strand=ignore.strand)
        coding[subjectHits(olap) %in% queryHits(hits)] <- TRUE
    } else {
        coding <- rep.int(NA, length(olap))
    }
    cmptrans <- .compatibleTranscription(query, subject, olap)
    compatible <- cmptrans$compatible
    unique <- .oneMatch(compatible, queryHits(olap))
    strandSpecific <- all(strand(query) != "*")[queryHits(olap)]
    if (!any(nc <- !compatible)) {
        ## compatible only 
        warning("all ranges in 'query' are compatible with ranges in 'subject'")
        .result(olap, compatible, unique, coding, strandSpecific)
    } else {
        ## compatible and non-compatible
        splice <- .gaps(query)
        intron <- .gaps(subject)
        novelExon <- .novelExon(splice, intron, nc)
        novelBounds <- .novelBounds(query, subject, olap)
        novelRetention <- .novelRetention(query, intron, nc)
        novelSpliceEvent <- .novelSpliceEvent(splice, intron)

        .result(olap, compatible, unique, coding, strandSpecific, nc=nc,
                novelTSS=novelBounds$TSS, novelTSE=novelBounds$TSE,
                novelSite=novelSpliceEvent$Site,
                novelJunction=novelSpliceEvent$Junction,
                novelExon, novelRetention)
    }
}

## -------------------------------------------------------------------------
## Methods in Rsamtools 
##

## findSpliceOverlaps,character,ANY-method
## findSpliceOverlaps,BamFile,ANY-method
 
## -------------------------------------------------------------------------
## Methods in GenomicFeatures 
##

## findSpliceOverlaps,GappedAlignments,TranscriptDb-method
## findSpliceOverlaps,GappedAlignmentPairs,TranscriptDb-method
## findSpliceOverlaps,GRangesList,TranscriptDb-method
