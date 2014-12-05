### =========================================================================
### findOverlaps methods
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### findOverlaps_GenomicRanges_old()
###

.putRangesOnFirstCircle <- function(x, circle.length)
{
    x_start0 <- start(x) - 1L  # 0-based start
    x_shift0 <- x_start0 %% circle.length - x_start0
    shift(x, x_shift0)
}

### 'circle.length' must be NA (if the underlying sequence is linear) or the
### length of the underlying circular sequence (integer vector of length 1
### with the name of the sequence).
### 'query' and 'subject' must be IRanges objects.
.findOverlaps.circle <- function(circle.length, query, subject,
                                 maxgap, minoverlap, type)
{
    if (is.na(circle.length))
        return(findOverlaps(query, subject,
                            maxgap=maxgap, minoverlap=minoverlap,
                            type=type, select="all",
                            algorithm="intervaltree"))
    q_len <- length(query)
    s_len <- length(subject)
    if (q_len == 0L || s_len == 0L)
        return(Hits(queryLength=q_len, subjectLength=s_len))
    query0 <- .putRangesOnFirstCircle(query, circle.length)
    subject0 <- .putRangesOnFirstCircle(subject, circle.length)
    pp_subject <- IntervalTree(subject0)
    hits00 <- findOverlaps(query0, pp_subject,
                           maxgap=maxgap, minoverlap=minoverlap,
                           type=type, select="all", algorithm="intervaltree")
    query1 <- shift(query0, circle.length)
    hits10 <- findOverlaps(query1, pp_subject,
                           maxgap=maxgap, minoverlap=minoverlap,
                           type=type, select="all", algorithm="intervaltree")
    subject1 <- shift(subject0, circle.length)
    hits01 <- findOverlaps(query0, subject1,
                           maxgap=maxgap, minoverlap=minoverlap,
                           type=type, select="all", algorithm="intervaltree")
    ## Merge 'hits00', 'hits10' and 'hits01'.
    union(union(hits00, hits10), hits01)
}

### 'x' must be 'strand(query)' or 'strand(subject)'.
.strandAsSignedNumber <- function(x)
{
    tmp <- as.integer(runValue(x))
    idx <- tmp >= 2L
    tmp[idx] <- tmp[idx] - 3L
    runValue(x) <- tmp
    as.vector(x)
}

findOverlaps_GenomicRanges_old <- function(query, subject,
             maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
{
    if (!isSingleNumber(maxgap) || maxgap < 0L)
        stop("'maxgap' must be a non-negative integer")
    type <- match.arg(type)
    select <- match.arg(select)

    ## merge() also checks that 'query' and 'subject' are based on the
    ## same reference genome.
    seqinfo <- merge(seqinfo(query), seqinfo(subject))

    q_len <- length(query)
    s_len <- length(subject)
    q_seqnames <- seqnames(query)
    s_seqnames <- seqnames(subject)
    q_splitranges <- splitRanges(q_seqnames)
    s_splitranges <- splitRanges(s_seqnames)
    q_seqlevels_nonempty <- names(q_splitranges)[sapply(q_splitranges, length) > 0]
    s_seqlevels_nonempty <- names(s_splitranges)[sapply(s_splitranges, length) > 0]
    q_ranges <- unname(ranges(query))
    s_ranges <- unname(ranges(subject))
    if (ignore.strand) {
        q_strand <- rep.int(1L, q_len)
        s_strand <- rep.int(1L, s_len)
    } else {
        q_strand <- .strandAsSignedNumber(strand(query))
        s_strand <- .strandAsSignedNumber(strand(subject))
    }

    common_seqlevels <- intersect(q_seqlevels_nonempty, s_seqlevels_nonempty)
    results <- lapply(common_seqlevels,
        function(seqlevel)
        {
            if (isCircular(seqinfo)[seqlevel] %in% TRUE) {
                circle.length <- seqlengths(seqinfo)[seqlevel]
            } else {
                circle.length <- NA
            }
            q_idx <- q_splitranges[[seqlevel]]
            s_idx <- s_splitranges[[seqlevel]]
            hits <- .findOverlaps.circle(circle.length,
                        extractROWS(q_ranges, q_idx),
                        extractROWS(s_ranges, s_idx),
                        maxgap, minoverlap, type)
            q_hits <- queryHits(hits)
            s_hits <- subjectHits(hits)
            compatible_strand <-
                extractROWS(q_strand, q_idx)[q_hits] *
                extractROWS(s_strand, s_idx)[s_hits] != -1L
            hits <- hits[compatible_strand]
            remapHits(hits, query.map=as.integer(q_idx),
                            new.queryLength=q_len,
                            subject.map=as.integer(s_idx),
                            new.subjectLength=s_len)
        })

    ## Combine the results.
    q_hits <- unlist(lapply(results, queryHits))
    if (is.null(q_hits))
        q_hits <- integer(0)

    s_hits <- unlist(lapply(results, subjectHits))
    if (is.null(s_hits))
        s_hits <- integer(0)

    selectHits(Hits(q_hits, s_hits, q_len, s_len), select=select)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "findOverlaps" methods for GenomicRanges objects
###

findOverlaps_GenomicRanges <- function(query, subject,
             maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
{
    min.score <- IRanges:::min_overlap_score(maxgap, minoverlap)
    type <- match.arg(type)
    select <- match.arg(select)
    algorithm <- match.arg(algorithm)
    if (algorithm != "nclist")
        warning("'algorithm' is ignored when 'query' or 'subject' ",
                "is a GNCList object")
    findOverlaps_GNCList(query, subject, min.score=min.score,
                         type=type, select=select,
                         ignore.strand=ignore.strand)
}

setMethod("findOverlaps", c("GNCList", "GenomicRanges"),
    findOverlaps_GenomicRanges
)

setMethod("findOverlaps", c("GenomicRanges", "GNCList"),
    findOverlaps_GenomicRanges
)

setMethod("findOverlaps", c("GenomicRanges", "GenomicRanges"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        type <- match.arg(type)
        select <- match.arg(select)
        algorithm <- match.arg(algorithm)
        if (algorithm == "nclist")
            FUN <- findOverlaps_GenomicRanges
        else
            FUN <- findOverlaps_GenomicRanges_old
        FUN(query, subject,
            maxgap=maxgap, minoverlap=minoverlap,
            type=type, select=select,
            ignore.strand=ignore.strand)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "findOverlaps" methods for GRangesList objects
###

.overlap_score <- function(hits, query, subject)
{
    q_ranges <- ranges(query)[queryHits(hits)]
    s_ranges <- ranges(subject)[subjectHits(hits)]
    1L + pmin(end(q_ranges), end(s_ranges)) -
         pmax(start(q_ranges), start(s_ranges))
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
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        if (!isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)
        algorithm <- match.arg(algorithm)

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
                              algorithm=algorithm,
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
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        if (!isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)
        algorithm <- match.arg(algorithm)

        unlisted_query <- unlist(query, use.names=FALSE)
        query_groups <- togroup(query)

        ans00 <- findOverlaps(unlisted_query, subject,
                              maxgap=maxgap,
                              type=type, select="all",
                              algorithm=algorithm,
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
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        if (!isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)
        algorithm <- match.arg(algorithm)

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
                              algorithm=algorithm,
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
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        findOverlaps(as(query, "GRanges"), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     algorithm=match.arg(algorithm),
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
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        findOverlaps(as(query, "GRanges"), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     algorithm=match.arg(algorithm),
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
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        findOverlaps(query, as(subject, "GRanges"),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     algorithm=match.arg(algorithm),
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
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        findOverlaps(query, as(subject, "GRanges"),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     algorithm=match.arg(algorithm),
                     ignore.strand=ignore.strand)
    }
)

setMethod("findOverlaps", c("RangedData", "GenomicRanges"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        ## Calls "findOverlaps" method for c("RangesList", "GenomicRanges")
        ## defined above.
        findOverlaps(ranges(query), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     algorithm=match.arg(algorithm),
                     ignore.strand=ignore.strand)
    }
)

setMethod("findOverlaps", c("RangedData", "GRangesList"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        ## Calls "findOverlaps" method for c("RangesList", "GRangesList")
        ## defined above.
        findOverlaps(ranges(query), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     algorithm=match.arg(algorithm),
                     ignore.strand=ignore.strand)
    }
)

setMethod("findOverlaps", c("GenomicRanges", "RangedData"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        ## Calls "findOverlaps" method for c("GenomicRanges", "RangesList")
        ## defined above.
        findOverlaps(query, ranges(subject),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     algorithm=match.arg(algorithm),
                     ignore.strand=ignore.strand)
    }
)

setMethod("findOverlaps", c("GRangesList", "RangedData"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        ## Calls "findOverlaps" method for c("GRangesList", "RangesList")
        ## defined above.
        findOverlaps(query, ranges(subject),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     algorithm=match.arg(algorithm),
                     ignore.strand=ignore.strand)
    }
)

### WARNING: Unlike most findOverlaps() methods, the methods for
### SummarizedExperiment below return a Hits object 'ans' that is *not*
### consistent with 'query' (or 'subject'), in the sense that 'queryHits(ans)'
### (or 'subjectHits(ans)') is not a valid index into 'query' (or 'subject')
### when 'query' (or 'subject') is a SummarizedExperiment object.

setMethod("findOverlaps", c("SummarizedExperiment", "Vector"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        findOverlaps(rowData(query), subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     algorithm=match.arg(algorithm),
                     ignore.strand=ignore.strand)
    }
)

setMethod("findOverlaps", c("Vector", "SummarizedExperiment"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        findOverlaps(query, rowData(subject),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     algorithm=match.arg(algorithm),
                     ignore.strand=ignore.strand)
    }
)

setMethod("findOverlaps", c("SummarizedExperiment", "SummarizedExperiment"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        findOverlaps(rowData(query), rowData(subject),
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type), select=match.arg(select),
                     algorithm=match.arg(algorithm),
                     ignore.strand=ignore.strand)
    }
)


### =========================================================================
### findOverlaps-based methods
### -------------------------------------------------------------------------

countOverlaps_GenomicRanges <- function(query, subject,
              maxgap=0L, minoverlap=1L,
              type=c("any", "start", "end", "within", "equal"),
              algorithm=c("nclist", "intervaltree"),
              ignore.strand=FALSE)
{
    min.score <- IRanges:::min_overlap_score(maxgap, minoverlap)
    type <- match.arg(type)
    algorithm <- match.arg(algorithm)
    if (algorithm != "nclist")
        warning("'algorithm' is ignored when 'query' or 'subject' ",
                "is a GNCList object")
    ans <- findOverlaps_GNCList(query, subject, min.score=min.score,
                                type=type, select="count",
                                ignore.strand=ignore.strand)
    names(ans) <- names(query)
    ans
}

setMethod("countOverlaps", c("GNCList", "GenomicRanges"),
    countOverlaps_GenomicRanges
)

setMethod("countOverlaps", c("GenomicRanges", "GNCList"),
    countOverlaps_GenomicRanges
)

countOverlaps.definition <- function(query, subject,
        maxgap=0L, minoverlap=1L,
        type=c("any", "start", "end", "within", "equal"),
        algorithm=c("nclist", "intervaltree"),
        ignore.strand=FALSE)
{
    counts <- queryHits(findOverlaps(query, subject, maxgap=maxgap,
                                     minoverlap=minoverlap,
                                     type=match.arg(type),
                                     algorithm=algorithm,
                                     ignore.strand=ignore.strand))
    structure(tabulate(counts, NROW(query)), names=names(query))
}

setMethod("countOverlaps", c("GenomicRanges", "GenomicRanges"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end", "within", "equal"),
             algorithm=c("nclist", "intervaltree"),
             ignore.strand=FALSE)
    {
        type <- match.arg(type)
        algorithm <- match.arg(algorithm)
        if (algorithm == "nclist") {
            countOverlaps_GenomicRanges(query, subject,
                                        maxgap=maxgap, minoverlap=minoverlap,
                                        type=type,
                                        ignore.strand=ignore.strand)
        } else {
            countOverlaps.definition(query, subject,
                                     maxgap=maxgap, minoverlap=minoverlap,
                                     type=type,
                                     algorithm=algorithm,
                                     ignore.strand=ignore.strand)
        }
    }
)

.signatures1 <- list(
    c("GenomicRanges", "Vector"),
    c("Vector", "GenomicRanges"),

    c("GRangesList", "Vector"),
    c("Vector", "GRangesList"),
    c("GRangesList", "GRangesList"),
    c("GRanges", "GRangesList"),
    c("GRangesList", "GRanges"),

    c("SummarizedExperiment", "Vector"),
    c("Vector", "SummarizedExperiment"),
    c("SummarizedExperiment", "SummarizedExperiment")
)

setMethods("countOverlaps", .signatures1, countOverlaps.definition)

overlapsAny.definition <- function(query, subject,
        maxgap=0L, minoverlap=1L,
        type=c("any", "start", "end", "within", "equal"),
        algorithm=c("nclist", "intervaltree"),
        ignore.strand=FALSE)
{
    !is.na(findOverlaps(query, subject,
                        maxgap=maxgap, minoverlap=minoverlap,
                        type=match.arg(type), select="arbitrary",
                        algorithm=match.arg(algorithm),
                        ignore.strand=ignore.strand))
}

subsetByOverlaps.definition1 <- function(query, subject,
        maxgap=0L, minoverlap=1L,
        type=c("any", "start", "end", "within", "equal"),
        algorithm=c("nclist", "intervaltree"),
        ignore.strand=FALSE)
{
    i <- overlapsAny(query, subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type),
                     algorithm=match.arg(algorithm),
                     ignore.strand=ignore.strand)
    extractROWS(query, i)
}

.subsetByOverlaps.definition2 <- function(query, subject,
        maxgap=0L, minoverlap=1L,
        type=c("any", "start", "end", "within", "equal"),
        algorithm=c("nclist", "intervaltree"),
        ignore.strand=FALSE)
{
    i <- overlapsAny(query, subject,
                     maxgap=maxgap, minoverlap=minoverlap,
                     type=match.arg(type),
                     algorithm=match.arg(algorithm),
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
    c("GRangesList", "RangedData"),
    c("SummarizedExperiment", "Vector"),
    c("Vector", "SummarizedExperiment"),
    c("SummarizedExperiment", "SummarizedExperiment")
)

for (sig in .signatures2) {
    setMethod("overlapsAny", sig, overlapsAny.definition)
    if (sig[1L] == "RangesList")
        setMethod("subsetByOverlaps", sig, .subsetByOverlaps.definition2)
    else
        setMethod("subsetByOverlaps", sig, subsetByOverlaps.definition1)
}

