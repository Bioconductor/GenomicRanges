### =========================================================================
### findOverlaps methods
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "findOverlaps" methods for GenomicRanges objects
###

findOverlaps_GenomicRanges <- function(query, subject,
             maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
{
    type <- match.arg(type)
    select <- match.arg(select)
    findOverlaps_GNCList(query, subject,
                         maxgap=maxgap, minoverlap=minoverlap,
                         type=type, select=select,
                         ignore.strand=ignore.strand)
}

setMethod("findOverlaps", c("GenomicRanges", "GenomicRanges"),
    findOverlaps_GenomicRanges
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "findOverlaps" methods for GRangesList objects
###

.overlap_score <- function(hits, query, subject)
{
    q_ranges <- ranges(query)[queryHits(hits)]
    s_ranges <- ranges(subject)[subjectHits(hits)]
    1L + pmin.int(end(q_ranges), end(s_ranges)) -
         pmax.int(start(q_ranges), start(s_ranges))
}

.aggregated_sum <- function(x, f1, f2)
{
    sm <- S4Vectors:::selfmatchIntegerPairs(f1, f2)
    S4Vectors:::tabulate2(sm, length(sm), weight=x)[sm]
}

setMethod("findOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        if (!isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)

        unlisted_query <- unlist(query, use.names=FALSE)
        query_groups <- togroup(query)
        unlisted_subject <- unlist(subject, use.names=FALSE)
        subject_groups <- togroup(subject)
 
        if (type == "start") {
            keep <- which(S4Vectors:::diffWithInitialZero(subject_groups) != 0L)
            unlisted_subject <- unlisted_subject[keep]
            subject_groups <- subject_groups[keep]
        } else if (type == "end") {
            keep <- end(subject@partitioning)[elementLengths(subject) > 0L]
            unlisted_subject <- unlisted_subject[keep]
            subject_groups <- subject_groups[keep]
        }
 
        ans00 <- findOverlaps(unlisted_query, unlisted_subject,
                              maxgap=maxgap,
                              type=type, select="all",
                              ignore.strand=ignore.strand)

        if (minoverlap > 1L) {
            score <- .overlap_score(ans00, unlisted_query, unlisted_subject)
            score <- .aggregated_sum(score, query_groups[queryHits(ans00)],
                                            subject_groups[subjectHits(ans00)])
            mcols(ans00) <- DataFrame(score=score)
        } 
        if (type == "within") {
            ans01 <- remapHits(ans00, subject.map=subject_groups,
                                      new.subjectLength=length(subject))
            ans11 <- remapHits(ans01, query.map=query_groups,
                                      new.queryLength=length(query),
                                      with.counts=TRUE)
            keep_idx <- which(mcols(ans11)[ , "counts"] ==
                              elementLengths(query)[queryHits(ans11)])
            mcols(ans11) <- NULL
            ans <- ans11[keep_idx]
        } else {
            ans <- remapHits(ans00, query.map=query_groups,
                                    new.queryLength=length(query),
                                    subject.map=subject_groups,
                                    new.subjectLength=length(subject))
        }
        if (minoverlap > 1L) {
            keep_idx <- which(mcols(ans)[ , "score"] >= minoverlap)
            mcols(ans) <- NULL
            ans <- ans[keep_idx]
        }
        selectHits(ans, select=select)
    }
)

setMethod("findOverlaps", c("GRangesList", "GenomicRanges"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        if (!isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)

        unlisted_query <- unlist(query, use.names=FALSE)
        query_groups <- togroup(query)

        ans00 <- findOverlaps(unlisted_query, subject,
                              maxgap=maxgap,
                              type=type, select="all",
                              ignore.strand=ignore.strand)

        if (minoverlap > 1L) {
            score <- .overlap_score(ans00, unlisted_query, subject)
            score <- .aggregated_sum(score, query_groups[queryHits(ans00)],
                                            subjectHits(ans00))
            mcols(ans00) <- DataFrame(score=score)
        }
        if (type == "within") {
            ans10 <- remapHits(ans00, query.map=query_groups,
                                      new.queryLength=length(query),
                                      with.counts=TRUE)
            keep_idx <- which(mcols(ans10)[ , "counts"] ==
                              elementLengths(query)[queryHits(ans10)])
            mcols(ans10) <- NULL
            ans <- ans10[keep_idx]
        } else {
            ans <- remapHits(ans00, query.map=query_groups,
                                    new.queryLength=length(query))
        }
        if (minoverlap > 1L) {
            keep_idx <- which(mcols(ans)[ , "score"] >= minoverlap)
            mcols(ans) <- NULL
            ans <- ans[keep_idx]
        }
        selectHits(ans, select=select)
    }
)

setMethod("findOverlaps", c("GenomicRanges", "GRangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        if (!isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)

        unlisted_subject <- unlist(subject, use.names=FALSE)
        subject_groups <- togroup(subject)

        if (type == "start") {
            keep <- which(S4Vectors:::diffWithInitialZero(subject_groups) != 0L)
            unlisted_subject <-  unlisted_subject[keep]
            subject_groups <- subject_groups[keep]
        } else if (type == "end") {
            keep <- end(subject@partitioning)[elementLengths(subject) > 0L]
            unlisted_subject <-  unlisted_subject[keep]
            subject_groups <- subject_groups[keep]
        }

        ans00 <- findOverlaps(query, unlisted_subject,
                              maxgap=maxgap,
                              type=type, select="all",
                              ignore.strand=ignore.strand)

        if(minoverlap > 1L) {
            score <- .overlap_score(ans00, query, unlisted_subject)
            score <- .aggregated_sum(score, queryHits(ans00),
                                            subject_groups[subjectHits(ans00)])
            mcols(ans00) <- DataFrame(score=score)
        }
        ans <- remapHits(ans00, subject.map=subject_groups,
                                new.subjectLength=length(subject))
        if (minoverlap > 1L) {
            keep_idx <- which(mcols(ans)[ , "score"] >= minoverlap)
            mcols(ans) <- NULL
            ans <- ans[keep_idx]
        }
        selectHits(ans, select=select)
    }
)

### Needed by chipseq package.
setMethod("findOverlaps", c("RangedData", "GenomicRanges"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        ## Calls "findOverlaps" method for c("GenomicRanges", "GenomicRanges")
        ## defined above.
        findOverlaps(as(query, "GRanges"), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     ignore.strand=ignore.strand)
    }
)


### =========================================================================
### findOverlaps-based methods
### -------------------------------------------------------------------------

countOverlaps_GenomicRanges <- function(query, subject,
              maxgap=0L, minoverlap=1L,
              type=c("any", "start", "end", "within", "equal"),
              ignore.strand=FALSE)
{
    type <- match.arg(type)
    ans <- findOverlaps_GNCList(query, subject,
                                maxgap=maxgap, minoverlap=minoverlap,
                                type=type, select="count",
                                ignore.strand=ignore.strand)
    names(ans) <- names(query)
    ans
}

setMethod("countOverlaps", c("GenomicRanges", "GenomicRanges"),
    countOverlaps_GenomicRanges
)

