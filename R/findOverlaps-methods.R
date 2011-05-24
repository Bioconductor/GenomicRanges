### =========================================================================
### findOverlaps methods
### -------------------------------------------------------------------------

### 'circle.length' must be NA (if the underlying sequence is linear) or the
### length of the underlying circular sequence (integer vector of length 1
### with the name of the sequence).
### 'query' and 'subject' must be IRanges objects.
.findOverlaps.circle <- function(circle.length, query, subject,
                                 maxgap, minoverlap, type)
{
    if (is.na(circle.length))
        return(findOverlaps(query, subject,
                            maxgap = maxgap, minoverlap = minoverlap,
                            type = type, select = "all"))
    if (type != "any")
        stop("overlap type \"", type, "\" is not yet supported ",
             "for circular sequence ", names(circle.length))
    subject.shift0 <- (start(subject) - 1L) %% circle.length +
                      1L - start(subject)
    subject0 <- shift(subject, subject.shift0)
    inttree0 <- IntervalTree(subject0)
    query.shift0 <- (start(query) - 1L) %% circle.length +
                    1L - start(query)
    query0 <- shift(query, query.shift0)
    overlaps00 <- findOverlaps(query0, inttree0,
                               maxgap = maxgap, minoverlap = minoverlap,
                               type = type, select = "all")
    query1 <- shift(query0, circle.length)
    overlaps10 <- findOverlaps(query1, inttree0,
                               maxgap = maxgap, minoverlap = minoverlap,
                               type = type, select = "all")
    subject1 <- shift(subject0, circle.length)
    overlaps01 <- findOverlaps(query0, subject1,
                               maxgap = maxgap, minoverlap = minoverlap,
                               type = type, select = "all")
    ## Merge 'overlaps00', 'overlaps10' and 'overlaps01'.
    qHits <- c(queryHits(overlaps00),
               queryHits(overlaps10),
               queryHits(overlaps01))
    sHits <- c(subjectHits(overlaps00),
               subjectHits(overlaps10),
               subjectHits(overlaps01))
    matchDataFrame <- data.frame(query = qHits, subject = sHits)
    matchDataFrame <- matchDataFrame[!duplicated(matchDataFrame), , drop=FALSE]
    row.names(matchDataFrame) <- NULL
    new("RangesMatching", matchMatrix = as.matrix(matchDataFrame),
                          DIM = overlaps00@DIM)
}

setMethod("findOverlaps", c("GenomicRanges", "GenomicRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end", "within", "equal"),
             select = c("all", "first"), ignore.strand = FALSE)
    {
        #        if (!identical(minoverlap, 1L))
        #    warning("'minoverlap' argument is ignored")
        if (!IRanges:::isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)

        ## merge() also checks that 'query' and 'subject' are based on the
        ## same reference genome.
        seqinfo <- merge(seqinfo(query), seqinfo(subject))
        DIM <- c(length(query), length(subject))
        if (min(DIM) == 0L) {
            matchMatrix <-
              matrix(integer(), ncol = 2,
                     dimnames = list(NULL, c("query", "subject")))
        } else {
            querySeqnames <- seqnames(query)
            querySplitRanges <- splitRanges(querySeqnames)
            uniqueQuerySeqnames <- names(querySplitRanges)

            subjectSeqnames <- seqnames(subject)
            subjectSplitRanges <- splitRanges(subjectSeqnames)
            uniqueSubjectSeqnames <- names(subjectSplitRanges)

            commonSeqnames <-
              intersect(uniqueQuerySeqnames, uniqueSubjectSeqnames)
             
            if (ignore.strand) {
                queryStrand <- rep.int(1L, length(query))
                subjectStrand <- rep.int(1L, length(subject))
            } else {
                queryStrand <- strand(query)
                levels(queryStrand) <- c("1", "-1", "0")
                queryStrand@values <-
                as.integer(as.character(runValue(queryStrand)))
                queryStrand <- as.vector(queryStrand)
                
                subjectStrand <- strand(subject)
                levels(subjectStrand) <- c("1", "-1", "0")
                subjectStrand@values <-
                as.integer(as.character(runValue(subjectStrand)))
                subjectStrand <- as.vector(subjectStrand)
            }
            queryRanges <- unname(ranges(query))
            subjectRanges <- unname(ranges(subject))

            matchMatrix <-
              do.call(rbind,
                      lapply(commonSeqnames, function(seqnm)
                      {
                          if (isCircular(seqinfo)[seqnm] %in% TRUE)
                              circle.length <- seqlengths(seqinfo)[seqnm]
                          else
                              circle.length <- NA
                          qIdxs <- querySplitRanges[[seqnm]]
                          sIdxs <- subjectSplitRanges[[seqnm]]
                          overlaps <- .findOverlaps.circle(
                                          circle.length,
                                          seqselect(queryRanges, qIdxs),
                                          seqselect(subjectRanges, sIdxs),
                                          maxgap, minoverlap, type)
                          qHits <- queryHits(overlaps)
                          sHits <- subjectHits(overlaps)
                          matches <-
                            cbind(query = as.integer(qIdxs)[qHits],
                                  subject = as.integer(sIdxs)[sHits])
                          matches[which(seqselect(queryStrand, qIdxs)[qHits] *
                                        seqselect(subjectStrand, sIdxs)[sHits] != -1L), ,
                                  drop=FALSE]
                      }))
            if (is.null(matchMatrix)) {
                matchMatrix <-
                  matrix(integer(0), nrow = 0, ncol = 2,
                         dimnames = list(NULL, c("query", "subject")))
            }
            matchMatrix <-
              matchMatrix[IRanges:::orderTwoIntegers(matchMatrix[ , 1L, drop=TRUE],
                                                     matchMatrix[ , 2L, drop=TRUE]), ,
                          drop=FALSE]
        }
        if (select == "all") {
            new("RangesMatching", matchMatrix = matchMatrix, DIM = DIM)
        } else {
            IRanges:::.matchMatrixToVector(matchMatrix, length(query))
        }
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
### IRanges:::duplicatedTwoIntegers() and then sort them. Could this be
### faster?
.cleanMatchMatrix <- function(matchMatrix)
{
    if (nrow(matchMatrix) <= 1L)
        return(matchMatrix)
    ## First sort the rows.
    oo <- IRanges:::orderTwoIntegers(matchMatrix[ , 1L, drop=TRUE],
                                     matchMatrix[ , 2L, drop=TRUE])
    matchMatrix <- matchMatrix[oo, , drop=FALSE]
    ## Then keep the unique rows.
    keep <- IRanges:::runEndsOfIntegerPairs(matchMatrix[ , 1L, drop=TRUE],
                                            matchMatrix[ , 2L, drop=TRUE])
    matchMatrix[keep, , drop=FALSE]
}

.groupSums <- function(x, by)
{
    f <- paste(by[,1], by[,2], sep="|")
    f <- factor(f, levels=unique(f))
    cil <- IRanges:::newCompressedList("CompressedIntegerList",
                       unlistData=x,
                       splitFactor=f)
    v <- Views(cil@unlistData, cil@partitioning)
    viewSums(v)
}

.updateMatchMatrix <- function(matchMatrix, intrsct, minoverlap) {
    widthSum <- .groupSums(width(intrsct), matchMatrix)
    dups <- IRanges:::duplicatedTwoIntegers(matchMatrix[,"query"],
        matchMatrix[,"subject"])
    indx <- (widthSum >= minoverlap)
    matchMatrix <- matchMatrix[!dups,  ,drop = FALSE]           
    matchMatrix <- matchMatrix[indx,  ,drop = FALSE]  
}

.makeGRL2GRmatchMatrix <- function(mm00, qpartitioning,
                                   type.is.within)
{
    query0 <- unname(mm00[ , "query"])
    subject0 <- unname(mm00[ , "subject"])
    oo <- IRanges:::orderTwoIntegers(subject0, query0)
    mm00 <- mm00[oo, , drop=FALSE]

    query1 <- togroup(qpartitioning, unname(mm00[ , "query"]))
    subject0 <- unname(mm00[ , "subject"])
    runend <- IRanges:::runEndsOfIntegerPairs(query1, subject0)
    mm10 <- cbind(query=query1, subject=subject0)[runend, , drop=FALSE]

    if (type.is.within) {
        runlen <- IRanges:::diffWithInitialZero(runend)
        keep <- width(qpartitioning)[mm10[ , "query"]] == runlen
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
        queryGroups <- togroup(query@partitioning)
        ans <- findOverlaps(unlistQuery, subject,
                            maxgap = maxgap, type = type, select = "all",
                            ignore.strand = ignore.strand)
        mm00 <- ans@matchMatrix
        if (minoverlap > 1L && nrow(mm00) > 0L) {
            query1 <- queryGroups[queryHits(ans)]
            subject0 <- unname(mm00[ , "subject"])
            mm10 <- cbind(query=query1, subject=subject0)
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
            DIM <- c(length(query), length(subject))
            initialize(ans, matchMatrix = mm10, DIM = DIM)
        } else {
            IRanges:::.matchMatrixToVector(mm10, length(query))
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
        subjectGroups <- togroup(subject@partitioning)
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
        matchMatrix <- ans@matchMatrix
        if(minoverlap > 1L && nrow(ans@matchMatrix) > 0) {    
            matchMatrix[,"subject"] <-  subjectGroups[subjectHits(ans)]
            intrsct <- pintersect(ranges(query)[queryHits(ans)],
                        ranges(unlistSubject[subjectHits(ans)]))
            matchMatrix <- .updateMatchMatrix(matchMatrix, intrsct, minoverlap)
        } else{
            matchMatrix[, 2L] <- subjectGroups[matchMatrix[, 2L, drop=TRUE]]
            matchMatrix <- .cleanMatchMatrix(matchMatrix)
        }
        if (select == "all") {
            DIM <- c(length(query), length(subject))
            initialize(ans, matchMatrix = matchMatrix, DIM = DIM)
        } else {
            IRanges:::.matchMatrixToVector(matchMatrix, length(query))
        }
    }
)

.makeGRL2GRLmatchMatrix <- function(mm00, qpartitioning, spartitioning,
                                    type.is.within)
{
    query0 <- unname(mm00[ , "query"])
    subject1 <- togroup(spartitioning, unname(mm00[ , "subject"]))
    oo <- IRanges:::orderTwoIntegers(subject1, query0)
    mm01 <- unique(cbind(query=query0, subject=subject1)[oo, , drop=FALSE])

    query1 <- togroup(qpartitioning, unname(mm01[ , "query"]))
    subject1 <- unname(mm01[ , "subject"])
    runend <- IRanges:::runEndsOfIntegerPairs(query1, subject1)
    mm11 <- cbind(query=query1, subject=subject1)[runend, , drop=FALSE]

    if (type.is.within) {
        runlen <- IRanges:::diffWithInitialZero(runend)
        keep <- width(qpartitioning)[mm11[ , "query"]] == runlen
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
        subjectGroups <- togroup(subject@partitioning)
        unlistQuery <- unlist(query, use.names = FALSE)
        queryGroups <- togroup(query@partitioning)
       
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
        mm00 <- ans@matchMatrix
        if (minoverlap > 1L && nrow(mm00) > 0L) {
            query1 <- queryGroups[queryHits(ans)]
            subject1 <- subjectGroups[subjectHits(ans)]
            mm11 <- cbind(query=query1, subject=subject1)
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
                                            type == "within")
        }
        ## Only for sorting, rows are already unique.
        ## TODO: Optimize this (.cleanMatchMatrix is also extracting unique
        ## rows but this is not necessary since they are already unique).
        mm11 <- .cleanMatchMatrix(mm11)
        if (select == "all") {
            DIM <- c(length(query), length(subject))
            initialize(ans, matchMatrix = mm11, DIM = DIM)
        } else {
            IRanges:::.matchMatrixToVector(mm11, length(query))
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

setMethod("findOverlaps", c("GappedAlignments", "ANY"),
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

setMethod("findOverlaps", c("ANY", "GappedAlignments"),
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
###   Note: Method with signature "GappedAlignments#ANY" chosen for
###    function "findOverlaps", target signature
###    "GappedAlignments#GappedAlignments".
###    "ANY#GappedAlignments" would also be valid
setMethod("findOverlaps", c("GappedAlignments", "GappedAlignments"),
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


### =========================================================================
### findOverlaps-based methods
### -------------------------------------------------------------------------

.countOverlaps.default <- function(query, subject,
        maxgap = 0L, minoverlap = 1L,
        type = c("any", "start", "end", "within", "equal"),
        ignore.strand = FALSE)
{
    tabulate(queryHits(findOverlaps(query, subject, maxgap = maxgap,
                                    minoverlap = minoverlap,
                                    type = match.arg(type),
                                    ignore.strand = ignore.strand)),
             NROW(query))
}

.subsetByOverlaps.default <- function(query, subject,
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

.subsetByOverlaps2 <- function(query, subject,
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

.subsetByOverlaps3 <- function(query, subject,
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

.match.default <- function(x, table,
                           nomatch = NA_integer_, incomparables = NULL)
{
    if (length(nomatch) != 1L)
        stop("'nomatch' must be of length 1")
    ans <- findOverlaps(x, table, select = "first")
    if (!is.na(nomatch) && IRanges:::anyMissing(ans))
        ans[is.na(ans)] <- nomatch
    ans
}

.signatures <- list(
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
    c("GappedAlignments", "ANY"),
    c("ANY", "GappedAlignments"),
    c("GappedAlignments", "GappedAlignments")
)

for (sig in .signatures) {
    setMethod("countOverlaps", sig, .countOverlaps.default)
    if (sig[1L] == "RangesList")
        setMethod("subsetByOverlaps", sig, .subsetByOverlaps2)
    else if (sig[1L] == "RangedData")
        setMethod("subsetByOverlaps", sig, .subsetByOverlaps3)
    else
        setMethod("subsetByOverlaps", sig, .subsetByOverlaps.default)
    setMethod("match", sig, .match.default)
    setMethod("%in%", sig, function(x, table) !is.na(match(x, table)))
}

