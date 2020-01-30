### =========================================================================
### findOverlaps methods
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "findOverlaps" methods for GenomicRanges objects
###

findOverlaps_GenomicRanges <- function(query, subject,
             maxgap=-1L, minoverlap=0L,
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

.overlapwidth <- function(hits, query, subject)
{
    q_ranges <- ranges(query)[queryHits(hits)]
    s_ranges <- ranges(subject)[subjectHits(hits)]
    ## TODO: Replace the code below by a call to
    ## poverlapWidth(q_ranges, s_ranges) when it's available.
    score <- pmin.int(end(q_ranges), end(s_ranges)) -
                     pmax.int(start(q_ranges), start(s_ranges)) + 1L
    pmax.int(score, 0L)
}

.aggregated_sum <- function(x, f1, f2)
{
    sm <- selfmatchIntegerPairs(f1, f2)
    S4Vectors:::tabulate2(sm, length(sm), weight=x)[sm]
}

setMethod("findOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        if (!isSingleNumber(minoverlap) || minoverlap < 0L)
            stop("'minoverlap' must be a single non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)

        unlisted_query <- unlist(query, use.names=FALSE)
        query_groups <- togroup(PartitioningByWidth(query))
        unlisted_subject <- unlist(subject, use.names=FALSE)
        subject_groups <- togroup(PartitioningByWidth(subject))
 
        if (type == "start") {
            keep <- which(S4Vectors:::diffWithInitialZero(subject_groups) != 0L)
            unlisted_subject <- unlisted_subject[keep]
            subject_groups <- subject_groups[keep]
        } else if (type == "end") {
            keep <- end(subject@partitioning)[elementNROWS(subject) > 0L]
            unlisted_subject <- unlisted_subject[keep]
            subject_groups <- subject_groups[keep]
        }
 
        ans00 <- findOverlaps(unlisted_query, unlisted_subject,
                              maxgap=maxgap,
                              type=type, select="all",
                              ignore.strand=ignore.strand)

        if (minoverlap > 0L) {
            owidth <- .overlapwidth(ans00, unlisted_query, unlisted_subject)
            owidth <- .aggregated_sum(owidth,
                                      query_groups[queryHits(ans00)],
                                      subject_groups[subjectHits(ans00)])
            mcols(ans00) <- DataFrame(owidth=owidth)
        } 
        if (type == "within" || type == "equal") {
            ans01 <- remapHits(ans00, Rnodes.remapping=subject_groups,
                                      new.nRnode=length(subject))
            ans11 <- remapHits(ans01, Lnodes.remapping=query_groups,
                                      new.nLnode=length(query),
                                      with.counts=TRUE)
            keep <- mcols(ans11)[ , "counts"] ==
                elementNROWS(query)[queryHits(ans11)]
            if (type == "equal") {
                ans10 <- remapHits(ans00, Lnodes.remapping=query_groups,
                                   new.nLnode=length(query))
                ans11 <- remapHits(ans10, Rnodes.remapping=subject_groups,
                                   new.nRnode=length(subject),
                                   with.counts=TRUE)
                keep <- keep & mcols(ans11)[ , "counts"] ==
                    elementNROWS(subject)[subjectHits(ans11)]
            }
            mcols(ans11) <- NULL
            ans <- ans11[keep]
        } else {
            ans <- remapHits(ans00, Lnodes.remapping=query_groups,
                                    new.nLnode=length(query),
                                    Rnodes.remapping=subject_groups,
                                    new.nRnode=length(subject))
        }
        if (minoverlap > 0L) {
            keep_idx <- which(mcols(ans, use.names=FALSE)[ , "owidth"] >=
                              minoverlap)
            mcols(ans) <- NULL
            ans <- ans[keep_idx]
        }
        selectHits(ans, select=select)
    }
)

setMethod("findOverlaps", c("GRangesList", "GenomicRanges"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        if (!isSingleNumber(minoverlap) || minoverlap < 0L)
            stop("'minoverlap' must be a single non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)

        unlisted_query <- unlist(query, use.names=FALSE)
        query_groups <- togroup(PartitioningByWidth(query))

        ans00 <- findOverlaps(unlisted_query, subject,
                              maxgap=maxgap,
                              type=type, select="all",
                              ignore.strand=ignore.strand)

        if (minoverlap > 0L) {
            owidth <- .overlapwidth(ans00, unlisted_query, subject)
            owidth <- .aggregated_sum(owidth,
                                      query_groups[queryHits(ans00)],
                                      subjectHits(ans00))
            mcols(ans00) <- DataFrame(owidth=owidth)
        }
        if (type == "within" || type == "equal") {
            ans10 <- remapHits(ans00, Lnodes.remapping=query_groups,
                                      new.nLnode=length(query),
                                      with.counts=TRUE)
            keep_idx <- which(mcols(ans10, use.names=FALSE)[ , "counts"] ==
                              elementNROWS(query)[queryHits(ans10)])
            mcols(ans10) <- NULL
            ans <- ans10[keep_idx]
        } else {
            ans <- remapHits(ans00, Lnodes.remapping=query_groups,
                                    new.nLnode=length(query))
        }
        if (minoverlap > 0L) {
            keep_idx <- which(mcols(ans, use.names=FALSE)[ , "owidth"] >=
                              minoverlap)
            mcols(ans) <- NULL
            ans <- ans[keep_idx]
        }
        selectHits(ans, select=select)
    }
)

setMethod("findOverlaps", c("GenomicRanges", "GRangesList"),
    function(query, subject, maxgap=-1L, minoverlap=0L,
             type=c("any", "start", "end", "within", "equal"),
             select=c("all", "first", "last", "arbitrary"),
             ignore.strand=FALSE)
    {
        if (!isSingleNumber(minoverlap) || minoverlap < 0L)
            stop("'minoverlap' must be a single non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)

        unlisted_subject <- unlist(subject, use.names=FALSE)
        subject_groups <- togroup(PartitioningByWidth(subject))

        if (type == "start") {
            keep <- which(S4Vectors:::diffWithInitialZero(subject_groups) != 0L)
            unlisted_subject <-  unlisted_subject[keep]
            subject_groups <- subject_groups[keep]
        } else if (type == "end") {
            keep <- end(subject@partitioning)[elementNROWS(subject) > 0L]
            unlisted_subject <-  unlisted_subject[keep]
            subject_groups <- subject_groups[keep]
        }

        ans00 <- findOverlaps(query, unlisted_subject,
                              maxgap=maxgap,
                              type=type, select="all",
                              ignore.strand=ignore.strand)

        if(minoverlap > 0L) {
            owidth <- .overlapwidth(ans00, query, unlisted_subject)
            owidth <- .aggregated_sum(owidth,
                                      queryHits(ans00),
                                      subject_groups[subjectHits(ans00)])
            mcols(ans00) <- DataFrame(owidth=owidth)
        }
        ans <- remapHits(ans00, Rnodes.remapping=subject_groups,
                                new.nRnode=length(subject),
                                with.counts=(type == "equal"))
        if (type == "equal") {
            keep <- mcols(ans)[ , "counts"] ==
                elementNROWS(subject)[subjectHits(ans)]
            mcols(ans) <- NULL
            ans <- ans[keep]
        }
        if (minoverlap > 0L) {
            keep_idx <- which(mcols(ans, use.names=FALSE)[ , "owidth"] >=
                              minoverlap)
            mcols(ans) <- NULL
            ans <- ans[keep_idx]
        }
        selectHits(ans, select=select)
    }
)


### =========================================================================
### findOverlaps-based methods
### -------------------------------------------------------------------------

countOverlaps_GenomicRanges <- function(query, subject,
              maxgap=-1L, minoverlap=0L,
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

compatibleStrand <- function(a, b) {
    a == "*" | b == "*" | a == b
}

setMethod("poverlaps", c("GenomicRanges", "GenomicRanges"),
          function(query, subject, maxgap=0L, minoverlap=1L,
                   type=c("any", "start", "end", "within", "equal"),
                   ignore.strand=FALSE)
{
    seqnames(query) == seqnames(subject) &
        (if (ignore.strand) TRUE
         else compatibleStrand(strand(query), strand(subject))) &
        poverlaps(ranges(query), ranges(subject), maxgap, minoverlap, type)
})
