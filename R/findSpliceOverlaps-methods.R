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
    .findSpliceOverlaps(query, subject, ignore.strand, singleEnd=FALSE, cds=cds)
})

setMethod("findSpliceOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, ignore.strand=FALSE, ..., cds=NULL)
{
    .findSpliceOverlaps(query, subject, ignore.strand, cds=cds)
})

.findSpliceOverlaps <- function(query, subject, ignore.strand=FALSE,
                                singleEnd=TRUE, cds=NULL)
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
    intronRegion <- .intronicRegions(subject, intron) 
    query <- query[queryHits(olap)]
    subject <- subject[subjectHits(olap)]
    intron <- intron[subjectHits(olap)]
    if (!singleEnd) {
        splice <- introns(query)
        query <- grglist(query, order.as.in.query=TRUE)
    } else { 
        splice <- .gaps(query)
    }
    compatible <- .compatibleTranscription(query, subject, splice)
    unique <- .oneMatch(compatible, queryHits(olap))
    strandSpecific <- all(strand(query) != "*")[queryHits(olap)]
    if (!any(nc <- !compatible)) {
        ## compatible only 
        warning("all ranges in 'query' are compatible with ranges in 'subject'")
        .result(olap, nc=NULL, compatible, unique, coding, strandSpecific)
    } else {
        ## non-compatible 
        novelTSS <- novelTSE <- novelExon <- novelRetention <- 
        novelSite <- novelJunction <- logical(length(olap))
        query <- query[nc]
        subject <- subject[nc]
        ## reads with splices
        novelExon[nc] <- .novelExon(splice[nc], intronRegion)
        novelSpliceEvent <- .novelSpliceEvent(splice[nc], intron[nc])
        novelSite[nc] <- novelSpliceEvent$Site 
        novelJunction[nc] <- novelSpliceEvent$Junction
        ## reads with or without splices
        novelBounds <- .novelBounds(query, subject, queryHits(olap)[nc])
        novelTSS[nc] <- novelBounds$TSS 
        novelTSE[nc] <- novelBounds$TSE
        novelRetention[nc] <- .novelRetention(query, intronRegion)
        .result(olap, compatible, unique, coding, strandSpecific, nc=nc,
                novelTSS, novelTSE, novelSite, novelJunction,
                novelExon, novelRetention)
    }
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
