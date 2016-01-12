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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other "findOverlaps" methods
###

### WARNING: Unlike most findOverlaps() methods, this method returns a Hits
### object 'ans' that is *not* consistent with 'query', in the sense that
### 'queryHits(ans)' is not a valid index into 'query'.
### Seems that the only use case for this method was to support the method
### for c("RangedData", "GenomicRanges").
### TODO: Deprecate after RangedData is gone.
setMethod("findOverlaps", c("RangesList", "GenomicRanges"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        findOverlaps(as(query, "GRanges"), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     ignore.strand=ignore.strand)
    }
)

### WARNING: Unlike most findOverlaps() methods, this method returns a Hits
### object 'ans' that is *not* consistent with 'query', in the sense that
### 'queryHits(ans)' is not a valid index into 'query'.
### Seems that the only use case for this method was to support the method
### for c("RangedData", "GRangesList").
### TODO: Deprecate after RangedData is gone.
setMethod("findOverlaps", c("RangesList", "GRangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        findOverlaps(as(query, "GRanges"), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     ignore.strand=ignore.strand)
    }
)

### WARNING: Unlike most findOverlaps() methods, this method returns a Hits
### object 'ans' that is *not* consistent with 'subject', in the sense that
### 'subjectHits(ans)' is not a valid index into 'subject'.
### Seems that the only use case for this method was to support the method
### for c("GenomicRanges", "RangedData").
### TODO: Deprecate after RangedData is gone.
setMethod("findOverlaps", c("GenomicRanges", "RangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        findOverlaps(query, as(subject, "GRanges"),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     ignore.strand=ignore.strand)
    }
)

### WARNING: Unlike most findOverlaps() methods, this method returns a Hits
### object 'ans' that is *not* consistent with 'subject', in the sense that
### 'subjectHits(ans)' is not a valid index into 'subject'.
### Seems that the only use case for this method was to support the method
### for c("GRangesList", "RangedData").
### TODO: Deprecate after RangedData is gone.
setMethod("findOverlaps", c("GRangesList", "RangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        findOverlaps(query, as(subject, "GRanges"),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     ignore.strand=ignore.strand)
    }
)

setMethod("findOverlaps", c("RangedData", "GenomicRanges"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        ## Calls "findOverlaps" method for c("RangesList", "GenomicRanges")
        ## defined above.
        findOverlaps(ranges(query), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     ignore.strand=ignore.strand)
    }
)

setMethod("findOverlaps", c("RangedData", "GRangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        ## Calls "findOverlaps" method for c("RangesList", "GRangesList")
        ## defined above.
        findOverlaps(ranges(query), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     ignore.strand=ignore.strand)
    }
)

setMethod("findOverlaps", c("GenomicRanges", "RangedData"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        ## Calls "findOverlaps" method for c("GenomicRanges", "RangesList")
        ## defined above.
        findOverlaps(query, ranges(subject),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     ignore.strand=ignore.strand)
    }
)

setMethod("findOverlaps", c("GRangesList", "RangedData"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        ## Calls "findOverlaps" method for c("GRangesList", "RangesList")
        ## defined above.
        findOverlaps(query, ranges(subject),
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

overlapsAny.definition <- function(query, subject,
        maxgap=0L, minoverlap=1L,
        type=c("any", "start", "end", "within", "equal"),
        ignore.strand=FALSE)
{
    !is.na(findOverlaps(query, subject,
                        maxgap=maxgap, minoverlap=minoverlap,
                        type=match.arg(type), select="arbitrary",
                        ignore.strand=ignore.strand))
}

subsetByOverlaps.definition1 <- function(query, subject,
        maxgap=0L, minoverlap=1L,
        type=c("any", "start", "end", "within", "equal"),
        ignore.strand=FALSE)
{
    i <- overlapsAny(query, subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type),
                     ignore.strand=ignore.strand)
    extractROWS(query, i)
}

.subsetByOverlaps.definition2 <- function(query, subject,
        maxgap=0L, minoverlap=1L,
        type=c("any", "start", "end", "within", "equal"),
        ignore.strand=FALSE)
{
    i <- overlapsAny(query, subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type),
                     ignore.strand=ignore.strand)
    query[splitAsList(i, space(query), drop=FALSE)]
}

.signatures2 <- list(
    c("GenomicRanges", "GenomicRanges"),
    c("GRangesList", "GenomicRanges"),
    c("GenomicRanges", "GRangesList"),
    c("GRangesList", "GRangesList"),
    c("RangesList", "GenomicRanges"),
    c("RangesList", "GRangesList"),
    c("GenomicRanges", "RangesList"),
    c("GRangesList", "RangesList"),
    c("RangedData", "GenomicRanges"),
    c("RangedData", "GRangesList"),
    c("GenomicRanges", "RangedData"),
    c("GRangesList", "RangedData")
)

for (sig in .signatures2) {
    setMethod("overlapsAny", sig, overlapsAny.definition)
    if (sig[1L] == "RangesList")
        setMethod("subsetByOverlaps", sig, .subsetByOverlaps.definition2)
    else
        setMethod("subsetByOverlaps", sig, subsetByOverlaps.definition1)
}

