### =========================================================================
### findSpliceOverlaps methods 
### -------------------------------------------------------------------------

setGeneric("findSpliceOverlaps", signature=c("query", "subject"),
    function(query, subject, ignore.strand=FALSE, ...)
{
    standardGeneric("findSpliceOverlaps")
})

### -------------------------------------------------------------------------
### subject as TranscriptDb 
###

#setMethod("findSpliceOverlaps", c("character", "TranscriptDb"),
#          function(query, subject, ignore.strand=FALSE, ...,
#                   param=ScanBamParam(), singleEnd=TRUE)
#{
#    callGeneric(BamFile(query), subject, ignore.strand, ...,
#                param=param, singleEnd=singleEnd)
#})
#
#setMethod("findSpliceOverlaps", c("BamFile", "TranscriptDb"),
#    function(query, subject, ignore.strand=FALSE, ...,
#             param=ScanBamParam(), singleEnd=TRUE)
#{
#    callGeneric(.readRanges(query, param, singleEnd), subject, 
#                ignore.strand, ...)
#})

setMethod("findSpliceOverlaps", c("GappedAlignments", "TranscriptDb"),
    function(query, subject, ignore.strand, ...)
{
    callGeneric(grglist(query, order.as.in.query=TRUE), subject, 
                ignore.strand, ...)
})

setMethod("findSpliceOverlaps", c("GappedAlignmentPairs", "TranscriptDb"),
    function(query, subject, ignore.strand, ...)
{
    callGeneric(grglist(query, order.as.in.query=TRUE), subject, 
                ignore.strand, ...)
})

setMethod("findSpliceOverlaps", c("GRangesList", "TranscriptDb"),
    function(query, subject, ignore.strand, ...)
{
    exbytx <- exonsBy(txdb, "tx")
    cds <- cdsBy(txdb, "tx")
    callGeneric(query, exbytx, ignore.strand, ..., cds=cds)
})

## -------------------------------------------------------------------------
## subject as Hits 
##

## -------------------------------------------------------------------------
## subject as GRangesList
##

#setMethod("findSpliceOverlaps", c("character", "GRangesList"),
#          function(query, subject, ignore.strand=FALSE, ..., 
#                   param=ScanBamParam(), singleEnd=TRUE, cds=NULL)
#{
#    callGeneric(BamFile(query), subject, ignore.strand, ..., 
#                param=param, singleEnd=singleEnd, cds=cds)
#})
#
#setMethod("findSpliceOverlaps", c("BamFile", "GRangesList"),
#    function(query, subject, ignore.strand=FALSE, ..., 
#             param=ScanBamParam(), singleEnd=TRUE, cds=NULL)
#{
#    callGeneric(.readRanges(query, param, singleEnd), subject, 
#                ignore.strand, ..., cds=cds)
#})

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
                novelSplicing=cmptrans$novelSplicing,
                novelTSS=novelBounds$TSS, novelTSE=novelBounds$TSE, 
                novelSite=novelSpliceEvent$Site, 
                novelJunction=novelSpliceEvent$Junction,
                novelExon, novelRetention)
    }
}

