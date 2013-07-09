### =========================================================================
### findOverlaps methods
### -------------------------------------------------------------------------

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
                            type=type, select="all"))
    q_len <- length(query)
    s_len <- length(subject)
    if (q_len == 0L || s_len == 0L)
        return(new("Hits", queryLength=q_len, subjectLength=s_len))
    if (type != "any")
        stop("overlap type \"", type, "\" is not yet supported ",
             "for circular sequence ", names(circle.length))
    subject0 <- .putRangesOnFirstCircle(subject, circle.length)
    inttree0 <- IntervalTree(subject0)
    query0 <- .putRangesOnFirstCircle(query, circle.length)
    hits00 <- findOverlaps(query0, inttree0,
                           maxgap=maxgap, minoverlap=minoverlap,
                           type=type, select="all")
    query1 <- shift(query0, circle.length)
    hits10 <- findOverlaps(query1, inttree0,
                           maxgap=maxgap, minoverlap=minoverlap,
                           type=type, select="all")
    subject1 <- shift(subject0, circle.length)
    hits01 <- findOverlaps(query0, subject1,
                           maxgap=maxgap, minoverlap=minoverlap,
                           type=type, select="all")
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

setMethod("findOverlaps", c("GenomicRanges", "GenomicRanges"),
    function(query, subject, maxgap=0L, minoverlap=1L,
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
        q_seqlevels <- levels(q_seqnames)
        s_seqlevels <- levels(s_seqnames)
        q_splitranges <- splitRanges(q_seqnames)
        s_splitranges <- splitRanges(s_seqnames)
        q_ranges <- unname(ranges(query))
        s_ranges <- unname(ranges(subject))
        if (ignore.strand) {
            q_strand <- rep.int(1L, q_len)
            s_strand <- rep.int(1L, s_len)
        } else {
            q_strand <- .strandAsSignedNumber(strand(query))
            s_strand <- .strandAsSignedNumber(strand(subject))
        }

        common_seqlevels <- intersect(q_seqlevels, s_seqlevels)
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
                                             seqselect(q_ranges, q_idx),
                                             seqselect(s_ranges, s_idx),
                                             maxgap, minoverlap, type)
                q_hits <- queryHits(hits)
                s_hits <- subjectHits(hits)
                compatible_strand <- seqselect(q_strand, q_idx)[q_hits] *
                                     seqselect(s_strand, s_idx)[s_hits] != -1L
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

        if (select == "arbitrary") {
            ans <- rep.int(NA_integer_, q_len)
            ans[q_hits] <- s_hits
            return(ans)
        }
        if (select == "first") {
            ans <- rep.int(NA_integer_, q_len)
            oo <- IRanges:::orderIntegerPairs(q_hits, s_hits, decreasing=TRUE)
            ans[q_hits[oo]] <- s_hits[oo]
            return(ans)
        }
        oo <- IRanges:::orderIntegerPairs(q_hits, s_hits)
        q_hits <- q_hits[oo]
        s_hits <- s_hits[oo]
        if (select == "last") {
            ans <- rep.int(NA_integer_, q_len)
            ans[q_hits] <- s_hits
            return(ans)
        }
        new2("Hits", queryHits=q_hits, subjectHits=s_hits,
                     queryLength=q_len, subjectLength=s_len,
                     check=FALSE)
    }
)

### Extract the unique rows from 2-col matrix 'matchMatrix' and return them
### sorted first by 1st col and then by 2nd col (actually the sorting is
### done first and then unique rows are extracted but this is an
### implementation detail since the final result doesn't depend on the
### order in which those things are done).
### TODO: Make this a function of 2 integer vectors (of equal length) and
### move it to IRanges/R/int-utils.R
### TODO: Try to invert the order i.e. first extract unique rows with
### IRanges:::duplicatedIntegerPairs() and then sort them. Could this be
### faster?
.cleanMatchMatrix <- function(matchMatrix)
{
    if (nrow(matchMatrix) <= 1L)
        return(matchMatrix)
    ## First sort the rows.
    oo <- IRanges:::orderIntegerPairs(matchMatrix[ , 1L],
                                      matchMatrix[ , 2L])
    matchMatrix <- matchMatrix[oo, , drop=FALSE]
    ## Then keep the unique rows.
    keep <- IRanges:::runEndsOfIntegerPairs(matchMatrix[ , 1L],
                                            matchMatrix[ , 2L])
    matchMatrix[keep, , drop=FALSE]
}

.groupSums <- function(x, by)
{
    f <- paste(by[,1], by[,2], sep="|")
    f <- factor(f, levels=unique(f))
    cil <- splitAsList(x, f)  # CompressedIntegerList
    v <- Views(cil@unlistData, cil@partitioning)
    viewSums(v)
}

.updateMatchMatrix <- function(matchMatrix, intrsct, minoverlap) {
    widthSum <- .groupSums(width(intrsct), matchMatrix)
    is_dup <- IRanges:::duplicatedIntegerPairs(matchMatrix[ , 1L],
                                               matchMatrix[ , 2L])
    indx <- (widthSum >= minoverlap)
    matchMatrix <- matchMatrix[!is_dup,  , drop=FALSE]           
    matchMatrix <- matchMatrix[indx,  , drop=FALSE]  
}

.makeGRL2GRmatchMatrix <- function(mm00, qpartitioning,
                                   type.is.within)
{
    query0 <- unname(mm00[ , 1L])
    subject0 <- unname(mm00[ , 2L])
    oo <- IRanges:::orderIntegerPairs(subject0, query0)
    mm00 <- mm00[oo, , drop=FALSE]

    query1 <- togroup(qpartitioning, j=unname(mm00[ , 1L]))
    subject0 <- unname(mm00[ , 2L])
    runend <- IRanges:::runEndsOfIntegerPairs(query1, subject0)
    mm10 <- cbind(queryHits=query1, subjectHits=subject0)[runend, , drop=FALSE]

    if (type.is.within) {
        runlen <- IRanges:::diffWithInitialZero(runend)
        keep <- width(qpartitioning)[mm10[ , 1L]] == runlen
        mm10 <- mm10[keep, , drop=FALSE]
    }
    mm10
}

setMethod("findOverlaps", c("GRangesList", "GenomicRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        if (!IRanges:::isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)
        unlistQuery <- unlist(query, use.names = FALSE)
        queryGroups <- togroup(query)
        ans <- findOverlaps(unlistQuery, subject,
                            maxgap = maxgap, type = type, select = "all",
                            ignore.strand = ignore.strand)
        mm00 <- as.matrix(ans)
        if (minoverlap > 1L && nrow(mm00) > 0L) {
            query1 <- queryGroups[queryHits(ans)]
            subject0 <- unname(mm00[ , 2L])
            mm10 <- cbind(queryHits=query1, subjectHits=subject0)
            intrsct <- pintersect(ranges(unlistQuery)[queryHits(ans)],
                                  ranges(subject)[subjectHits(ans)])
            mm10 <- .updateMatchMatrix(mm10, intrsct, minoverlap)
            if (type == "within") {
                ## TODO: Call .makeGRL2GRmatchMatrix() and intersect the
                ## result with 'mm10'.
                stop("'type=\"within\"' is not yet supported ",
                     "when 'minoverlap' > 1")
            }
        } else {
            mm10 <- .makeGRL2GRmatchMatrix(mm00,
                                           query@partitioning,
                                           type == "within")
        }
        ## Only for sorting, rows are already unique.
        ## TODO: Optimize this (.cleanMatchMatrix is also extracting unique
        ## rows but this is not necessary since they are already unique).
        mm10 <- .cleanMatchMatrix(mm10)
        if (select == "all") {
            initialize(ans,
                       queryHits = unname(mm10[ , 1L]),
                       subjectHits = unname(mm10[ , 2L]),
                       queryLength = length(query),
                       subjectLength = length(subject))
        } else {
            IRanges:::.hitsMatrixToVector(mm10, length(query))
        }
    }
)

setMethod("findOverlaps", c("GenomicRanges", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        if (!IRanges:::isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)

        unlistSubject <- unlist(subject, use.names=FALSE)
        subjectGroups <- togroup(subject)
        if (type == "start") {
            keep <- which(IRanges:::diffWithInitialZero(subjectGroups) != 0L)
            unlistSubject <-  unlistSubject[keep]
            subjectGroups <- subjectGroups[keep]
        } else if (type == "end") {
            keep <- end(subject@partitioning)[elementLengths(subject) > 0L]
            unlistSubject <-  unlistSubject[keep]
            subjectGroups <- subjectGroups[keep]
        }
        ans <- findOverlaps(query, unlistSubject,
                            maxgap = maxgap, type = type, select = "all",
                            ignore.strand = ignore.strand)
        matchMatrix <- as.matrix(ans)
        if(minoverlap > 1L && nrow(matchMatrix) > 0) {
            matchMatrix[ , 2L] <-  subjectGroups[subjectHits(ans)]
            intrsct <- pintersect(ranges(query)[queryHits(ans)],
                        ranges(unlistSubject[subjectHits(ans)]))
            matchMatrix <- .updateMatchMatrix(matchMatrix, intrsct, minoverlap)
        } else{
            matchMatrix[ , 2L] <- subjectGroups[matchMatrix[ , 2L]]
            matchMatrix <- .cleanMatchMatrix(matchMatrix)
        }
        if (select == "all") {
            initialize(ans,
                       queryHits = unname(matchMatrix[ , 1L]),
                       subjectHits = unname(matchMatrix[ , 2L]),
                       queryLength = length(query),
                       subjectLength = length(subject))
        } else {
            IRanges:::.hitsMatrixToVector(matchMatrix, length(query))
        }
    }
)

.makeGRL2GRLmatchMatrix <- function(mm00, qpartitioning, spartitioning,
                                    type.is.within)
{
    query0 <- unname(mm00[ , 1L])
    subject1 <- togroup(spartitioning, j=unname(mm00[ , 2L]))
    oo <- IRanges:::orderIntegerPairs(subject1, query0)
    mm01 <- cbind(queryHits=query0, subjectHits=subject1)[oo, , drop=FALSE]
    is_dup <- IRanges:::duplicatedIntegerPairs(mm01[ , 1L],
                                               mm01[ , 2L])
    mm01 <- mm01[!is_dup, , drop=FALSE]

    query1 <- togroup(qpartitioning, j=unname(mm01[ , 1L]))
    subject1 <- unname(mm01[ , 2L])
    runend <- IRanges:::runEndsOfIntegerPairs(query1, subject1)
    mm11 <- cbind(queryHits=query1, subjectHits=subject1)[runend, , drop=FALSE]

    if (type.is.within) {
        runlen <- IRanges:::diffWithInitialZero(runend)
        keep <- width(qpartitioning)[mm11[ , 1L]] == runlen
        mm11 <- mm11[keep, , drop=FALSE]
    }
    mm11
}

setMethod("findOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        if (!IRanges:::isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)

        unlistSubject <- unlist(subject, use.names=FALSE)
        subjectGroups <- togroup(subject)
        unlistQuery <- unlist(query, use.names = FALSE)
        queryGroups <- togroup(query)
 
        if (type == "start") {
            keep <- which(IRanges:::diffWithInitialZero(subjectGroups) != 0L)
            unlistSubject <- unlistSubject[keep]
            subjectGroups <- subjectGroups[keep]
        } else if (type == "end") {
            keep <- end(subject@partitioning)[elementLengths(subject) > 0L]
            unlistSubject <- unlistSubject[keep]
            subjectGroups <- subjectGroups[keep]
        }
 
        ans <- findOverlaps(unlistQuery, unlistSubject,
                            maxgap = maxgap, type = type, select = "all",
                            ignore.strand = ignore.strand)
        mm00 <- as.matrix(ans)
        if (minoverlap > 1L && nrow(mm00) > 0L) {
            query1 <- queryGroups[queryHits(ans)]
            subject1 <- subjectGroups[subjectHits(ans)]
            mm11 <- cbind(queryHits=query1, subjectHits=subject1)
            intrsct <- pintersect(ranges(unlistQuery)[queryHits(ans)],
                                  ranges(unlistSubject)[subjectHits(ans)])
            mm11 <- .updateMatchMatrix(mm11, intrsct, minoverlap)
            if (type == "within") {
                ## TODO: Call .makeGRL2GRLmatchMatrix() and intersect the
                ## result with 'mm11'.
                stop("'type=\"within\"' is not yet supported ",
                     "when 'minoverlap' > 1")
            }
        } else {
            mm11 <- .makeGRL2GRLmatchMatrix(mm00,
                                            query@partitioning,
                                            subject@partitioning,
                                            type.is.within = type == "within")
        }
        ## Only for sorting, rows are already unique.
        ## TODO: Optimize this (.cleanMatchMatrix is also extracting unique
        ## rows but this is not necessary since they are already unique).
        mm11 <- .cleanMatchMatrix(mm11)
        if (select == "all") {
            initialize(ans,
                       queryHits = unname(mm11[ , 1L]),
                       subjectHits = unname(mm11[ , 2L]),
                       queryLength = length(query),
                       subjectLength = length(subject))
        } else {
            IRanges:::.hitsMatrixToVector(mm11, length(query))
        }
    }
)

setMethod("findOverlaps", c("RangesList", "GenomicRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(as(query, "GRanges"), subject = subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("RangesList", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(as(query, "GRanges"), subject = subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GenomicRanges", "RangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(query, as(subject, "GRanges"),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GRangesList", "RangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(query, as(subject, "GRanges"),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("RangedData", "GenomicRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        ## Calls "findOverlaps" method for c("RangesList", "GenomicRanges")
        ## defined above.
        findOverlaps(ranges(query), subject = subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("RangedData", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        ## Calls "findOverlaps" method for c("RangesList", "GRangesList")
        ## defined above.
        findOverlaps(ranges(query), subject = subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GenomicRanges", "RangedData"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        ## Calls "findOverlaps" method for c("GenomicRanges", "RangesList")
        ## defined above.
        findOverlaps(query, ranges(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GRangesList", "RangedData"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        ## Calls "findOverlaps" method for c("GRangesList", "RangesList")
        ## defined above.
        findOverlaps(query, ranges(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GAlignments", "Vector"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(grglist(query), subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("Vector", "GAlignments"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(query, grglist(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

### Not strictly needed! Defining the above 2 methods covers that case but
### with the following note:
###   > findOverlaps(al1, al0)
###   Note: Method with signature "GAlignments#ANY" chosen for
###    function "findOverlaps", target signature
###    "GAlignments#GAlignments".
###    "ANY#GAlignments" would also be valid
setMethod("findOverlaps", c("GAlignments", "GAlignments"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(grglist(query), grglist(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GAlignments", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(grglist(query), subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GRangesList", "GAlignments"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(query, grglist(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GAlignmentPairs", "Vector"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(as(query, "GRangesList"), subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("Vector", "GAlignmentPairs"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(query, as(subject, "GRangesList"),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GAlignmentPairs", "GAlignmentPairs"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(as(query, "GRangesList"), as(subject, "GRangesList"),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("SummarizedExperiment", "Vector"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(rowData(query), subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("Vector", "SummarizedExperiment"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(query, rowData(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("SummarizedExperiment", "SummarizedExperiment"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        findOverlaps(rowData(query), rowData(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select),
                     ignore.strand = ignore.strand)
    }
)

setMethod("findOverlaps", c("GAlignmentsList", "Vector"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        hits <- findOverlaps(grglist(unlist(query, use.names = FALSE)),
                             subject, maxgap = maxgap, minoverlap = minoverlap,
                             type = match.arg(type), select = match.arg(select),
                             ignore.strand = ignore.strand)
        remapHits(hits, query.map=factor(togroup(query)))
    }
)

setMethod("findOverlaps", c("Vector", "GAlignmentsList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        hits <- findOverlaps(query, grglist(unlist(subject, use.names = FALSE)),
                             maxgap = maxgap, minoverlap = minoverlap,
                             type = match.arg(type), select = match.arg(select),
                             ignore.strand = ignore.strand)
        remapHits(hits, subject.map=factor(togroup(subject)))
    }
)

setMethod("findOverlaps", c("GAlignmentsList", "GAlignmentsList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        hits <- findOverlaps(grglist(unlist(query, use.names = FALSE)), 
                             grglist(unlist(subject, use.names = FALSE)),
                             maxgap = maxgap, minoverlap = minoverlap,
                             type = match.arg(type), 
                             select = match.arg(select),
                             ignore.strand = ignore.strand)
        remapHits(hits, subject.map=factor(togroup(subject)),
                  query.map=factor(togroup(query)))
    }
)

### =========================================================================
### findOverlaps-based methods
### -------------------------------------------------------------------------

.countOverlaps.definition <- function(query, subject,
        maxgap = 0L, minoverlap = 1L,
        type = c("any", "start", "end", "within", "equal"),
        ignore.strand = FALSE)
{
    counts <- queryHits(findOverlaps(query, subject, maxgap = maxgap,
                                     minoverlap = minoverlap,
                                     type = match.arg(type),
                                     ignore.strand = ignore.strand))
    structure(tabulate(counts, NROW(query)), names=names(query))
}

.signatures1 <- list(
    c("GenomicRanges", "Vector"),
    c("Vector", "GenomicRanges"),
    c("GenomicRanges", "GenomicRanges"),

    c("GRangesList", "Vector"),
    c("Vector", "GRangesList"),
    c("GRangesList", "GRangesList"),
    c("GRanges", "GRangesList"),
    c("GRangesList", "GRanges"),

    c("GAlignments", "Vector"),
    c("Vector", "GAlignments"),
    c("GAlignments", "GAlignments"),
    c("GAlignments", "GenomicRanges"),
    c("GenomicRanges", "GAlignments"),
    c("GAlignments", "GRangesList"),
    c("GRangesList", "GAlignments"),

    c("GAlignmentPairs", "Vector"),
    c("Vector", "GAlignmentPairs"),
    c("GAlignmentPairs", "GAlignmentPairs"),

    c("GAlignmentsList", "Vector"),
    c("Vector", "GAlignmentsList"),
    c("GAlignmentsList", "GAlignmentsList"),

    c("SummarizedExperiment", "Vector"),
    c("Vector", "SummarizedExperiment"),
    c("SummarizedExperiment", "SummarizedExperiment")
)

setMethods("countOverlaps", .signatures1, .countOverlaps.definition)

.overlapsAny.definition <- function(query, subject,
        maxgap = 0L, minoverlap = 1L,
        type = c("any", "start", "end", "within", "equal"),
        ignore.strand = FALSE)
{
    !is.na(findOverlaps(query, subject, maxgap = maxgap,
                        minoverlap = minoverlap,
                        type = match.arg(type),
                        select = "first",
                        ignore.strand = ignore.strand))
}

.subsetByOverlaps.definition1 <- function(query, subject,
        maxgap = 0L, minoverlap = 1L,
        type = c("any", "start", "end", "within", "equal"), 
        ignore.strand = FALSE)
{
    query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                              minoverlap = minoverlap,
                              type = match.arg(type),
                              select = "first",
                              ignore.strand = ignore.strand))]
}

.subsetByOverlaps.definition2 <- function(query, subject,
        maxgap = 0L, minoverlap = 1L,
        type = c("any", "start", "end", "within", "equal"), 
        ignore.strand = FALSE)
{
    i <- !is.na(findOverlaps(query, subject, maxgap = maxgap,
                             minoverlap = minoverlap,
                             type = match.arg(type),
                             select = "first",
                             ignore.strand = ignore.strand))
    query[seqsplit(i, space(query), drop=FALSE)]
}

.subsetByOverlaps.definition3 <- function(query, subject,
        maxgap = 0L, minoverlap = 1L,
        type = c("any", "start", "end", "within", "equal"), 
        ignore.strand = FALSE)
{
    query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                              minoverlap = minoverlap,
                              type = match.arg(type),
                              select = "first",
                              ignore.strand = ignore.strand)),]
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
    c("GAlignments", "Vector"),
    c("Vector", "GAlignments"),
    c("GAlignments", "GAlignments"),
    c("GAlignmentPairs", "Vector"),
    c("Vector", "GAlignmentPairs"),
    c("SummarizedExperiment", "Vector"),
    c("Vector", "SummarizedExperiment"),
    c("SummarizedExperiment", "SummarizedExperiment"),
    c("GAlignmentsList", "Vector"),
    c("Vector", "GAlignmentsList"),
    c("GAlignmentsList", "GAlignmentsList")
)

for (sig in .signatures2) {
    setMethod("overlapsAny", sig, .overlapsAny.definition)
    if (sig[1L] == "RangesList")
        setMethod("subsetByOverlaps", sig, .subsetByOverlaps.definition2)
    else if (sig[1L] %in% c("RangedData", "SummarizedExperiment"))
        setMethod("subsetByOverlaps", sig, .subsetByOverlaps.definition3)
    else
        setMethod("subsetByOverlaps", sig, .subsetByOverlaps.definition1)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### match(), %in%, findMatches(), and countMatches()
###
### All "match" and "%in% methods (except methods for GenomicRanges objects)
### are defunct in BioC 2.13.
### Methods for GenomicRanges have been moved to GenomicRanges-comparison.R
### and their behavior changed to do equality instead of overlaps.

.signatures3 <- .signatures2[-1L]

.match.definition <- function(x, table,
                              nomatch = NA_integer_, incomparables = NULL)
{
    if (!identical(nomatch, NA_integer_))
        stop("'nomatch' arg is not supported")
    msg <- c("match() between a ", class(x), " and a ", class(table),
             " object is defunct.\nPlease use '",
             "findOverlaps(x, table, select=\"first\")",
             "' instead.")
    .Defunct(msg=msg)  # deprecated in BioC 2.12, defunct in BioC 2.13
}

setMethods("match", .signatures3, .match.definition)

setMethods("%in%", .signatures3, IRanges:::`.%in%.definition`)

### The only reason for defining the methods below is to prevent the default
### "findMatches" or "countMatches" methods to be called and return something
### wrong (and the reason they would return something wrong is because they
### are based on match() which does overlaps instead of equality).
### TODO: Remove these methods in BioC 2.14 when the "match" methods for all
### the signatures in '.signatures3' are gone.

setMethods("findMatches", .signatures3,
    function(x, table, select=c("all", "first", "last"), ...)
    {
        msg <- c("findMatches() between a ", class(x), " and a ",
                 class(table), " object is not supported")
        stop(msg)
    }
)

setMethods("countMatches", .signatures3,
    function(x, table, ...)
    {
        msg <- c("countMatches() between a ", class(x), " and a ",
                 class(table), " object is not supported")
        stop(msg)
    }
)

