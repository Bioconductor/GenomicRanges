## -------------------------------------------------------------------------
## summarizeOverlaps methods and related utilities
## -------------------------------------------------------------------------

setGeneric("summarizeOverlaps", signature=c("features", "reads"),
    function(features, reads, mode=Union, ignore.strand=FALSE, ...)
{
    standardGeneric("summarizeOverlaps")
})

## -------------------------------------------------------------------------
## BamFiles and BamViews methods are in Rsamtools
##

## -------------------------------------------------------------------------
## GAlignments, GAlignmentsList and GAlignmentPairs methods
##

.dispatchOverlaps <-
    function(features, reads, mode, ignore.strand, inter.feature, ...)
{
    if (ignore.strand) {
        if (class(features) == "GRangesList") {
            r <- unlist(features)
            strand(r) <- "*"
            features@unlistData <- r
        } else {
            strand(features) <- "*"
        } 
    }
    mode(features, reads, ignore.strand, inter.feature=inter.feature)
}

.summarizeOverlaps <- function(features, reads, mode, ignore.strand,
                               ..., inter.feature=inter.feature)
{
    if (all(strand(reads) == "*"))
        ignore.strand <- TRUE
    mode <- match.fun(mode)
    counts <- .dispatchOverlaps(features, reads, mode, ignore.strand,
                                inter.feature=inter.feature)
    colData <- DataFrame(object=class(reads),
                         records=length(reads),
                         row.names="reads") 
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
                         rowData=features, colData=colData)
}

setMethod("summarizeOverlaps", c("GRanges", "GAlignments"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, reads, mode, ignore.strand,
                       ..., inter.feature=inter.feature)
})
setMethod("summarizeOverlaps", c("GRangesList", "GAlignments"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, reads, mode, ignore.strand,
                       ..., inter.feature=inter.feature)
})

setMethod("summarizeOverlaps", c("GRanges", "GAlignmentsList"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, grglist(reads, ignore.strand=TRUE), 
                       mode, ignore.strand, ..., 
                       inter.feature=inter.feature)
})
setMethod("summarizeOverlaps", c("GRangesList", "GAlignmentsList"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, grglist(reads, ignore.strand=TRUE), 
                       mode, ignore.strand, ..., 
                       inter.feature=inter.feature)
})

setMethod("summarizeOverlaps", c("GRanges", "GAlignmentPairs"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, grglist(reads), mode, ignore.strand,
                       ..., inter.feature=inter.feature)
})
setMethod("summarizeOverlaps", c("GRangesList", "GAlignmentPairs"),
    function(features, reads, mode, ignore.strand=FALSE, ...,
             inter.feature=TRUE)
{
    .summarizeOverlaps(features, grglist(reads), mode, ignore.strand,
                       ..., inter.feature=inter.feature)
})

## -------------------------------------------------------------------------
## 'mode' functions 
##

Union <- function(features, reads, ignore.strand=FALSE, inter.feature=TRUE)
{
    ov <- findOverlaps(features, reads, ignore.strand=ignore.strand)
    if (inter.feature) {
        ## Remove ambigous reads.
        reads_to_keep <- which(countSubjectHits(ov) == 1L)
        ov <- ov[subjectHits(ov) %in% reads_to_keep]
    }
    countQueryHits(ov)
}

## Drop from 'reads' circular seqlevels that are in use in *both*: 'reads'
## and 'features'.
.dropCircularSeqlevelsInUse <- function(reads, features)
{
    seqlevels_in_use <- intersect(seqlevelsInUse(reads),
                                  seqlevelsInUse(features))
    seqinfo <- merge(seqinfo(reads), seqinfo(features))
    is_circ <- isCircular(seqinfo)
    circular_seqlevels <- names(is_circ)[is_circ]
    seqlevels_to_drop <- intersect(seqlevels_in_use, circular_seqlevels)
    if (length(seqlevels_to_drop) != 0L) {
        warning("reads on circular sequence(s) '",
                paste(seqlevels_to_drop, sep="', '"),
                "' were ignored")
        seqlevels(reads, force=TRUE) <- setdiff(seqlevels(reads),
                                                seqlevels_to_drop)
    }
    reads
}

IntersectionStrict <- function(features, reads, ignore.strand=FALSE, 
                               inter.feature=TRUE)
{
    reads <- .dropCircularSeqlevelsInUse(reads, features)
    ov <- findOverlaps(reads, features, type="within",
                       ignore.strand=ignore.strand)
    if (inter.feature) {
        ## Remove ambigous reads.
        reads_to_keep <- which(countQueryHits(ov) == 1L)
        ov <- ov[queryHits(ov) %in% reads_to_keep]
    }
    countSubjectHits(ov)
}

.removeSharedRegions <- function(features, ignore.strand=FALSE)
{
    if (is(features, "GRanges")) {
        regions <- disjoin(features, ignore.strand=ignore.strand)
    } else if (is(features, "GRangesList")) {
        regions <- disjoin(features@unlistData, ignore.strand=ignore.strand)
    } else {
        stop("internal error")  # should never happen
    }
    ov <- findOverlaps(features, regions, ignore.strand=ignore.strand)
    regions_to_keep <- which(countSubjectHits(ov) == 1L)
    ov <- ov[subjectHits(ov) %in% regions_to_keep]
    ans_flesh <- regions[subjectHits(ov)]
    ### Using 'countQueryHits(ov)' to compute the skeleton of the answer relies
    ### on the assumption that the hits returned by findOverlaps() are always
    ### ordered by query first and then by subject.
    ans_eltlens <- countQueryHits(ov)
    ans_skeleton <- PartitioningByEnd(cumsum(ans_eltlens))
    relist(ans_flesh, ans_skeleton)
}

IntersectionNotEmpty <-  function(features, reads, ignore.strand=FALSE, 
                                  inter.feature=TRUE)
{
    features <- .removeSharedRegions(features, ignore.strand=ignore.strand)
    Union(features, reads, ignore.strand=ignore.strand,
           inter.feature=inter.feature)
}
