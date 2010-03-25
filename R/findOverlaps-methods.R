.similarSeqnameConvention <- function(seqs1, seqs2) {
    funList <-
      list(isRoman = function(x) grepl("[ivIV]", x),
           isArabic = function(x) grepl("[1-9]", x),
           haschr = function(x) grepl("^chr", x),
           hasChr = function(x) grepl("^Chr", x),
           hasCHR = function(x) grepl("^CHR", x),
           haschrom = function(x) grepl("^chrom", x),
           hasChrom = function(x) grepl("^Chrom", x),
           hasCHROM = function(x) grepl("^CHROM", x))
    all(sapply(funList, function(f) any(f(seqs1)) == any(f(seqs2))))
}


### =========================================================================
### findOverlaps methods
### -------------------------------------------------------------------------

.cleanMatchMatrix <- function(matchMatrix) {
    fastDiff <- IRanges:::diffWithInitialZero
    nr <- nrow(matchMatrix)
    nc <- ncol(matchMatrix)
    if (nr <= 1L) {
        matchMatrix
    } else {
        matchMatrix <-
          matchMatrix[IRanges:::orderTwoIntegers(matchMatrix[ , 1L, drop=TRUE],
                                                 matchMatrix[ , 2L, drop=TRUE]), ,
                      drop=FALSE]
        matchMatrix[fastDiff(matchMatrix[,1L,drop=TRUE]) != 0L |
                    fastDiff(matchMatrix[,2L,drop=TRUE]) != 0L, , drop=FALSE]
    }
}

.matchMatrixToVector <- function(matchMatrix, lengthQuery) {
    matchMatrix <-
      matchMatrix[IRanges:::diffWithInitialZero(matchMatrix[,1L,drop=TRUE]) != 0L,,
                  drop=FALSE]
    ans <- rep.int(NA_integer_, lengthQuery)
    ans[matchMatrix[,1L,drop=TRUE]] <- matchMatrix[,2L,drop=TRUE]
    ans
}

setMethod("findOverlaps", c("GRanges", "GRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        if (!IRanges:::isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)

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

            if (!.similarSeqnameConvention(uniqueQuerySeqnames,
                                           uniqueSubjectSeqnames))
                stop("'query' and 'subject' do not use a similiar naming ",
                     "convention for seqnames")

            commonSeqnames <-
              intersect(uniqueQuerySeqnames, uniqueSubjectSeqnames)

            queryStrand <- strand(query)
            levels(queryStrand) <- c("1", "-1", "0")
            queryStrand@values <-
              as.integer(as.character(runValue(queryStrand)))
            queryStrand <- as.vector(queryStrand)
            queryRanges <- unname(ranges(query))

            subjectStrand <- strand(subject)
            levels(subjectStrand) <- c("1", "-1", "0")
            subjectStrand@values <-
              as.integer(as.character(runValue(subjectStrand)))
            subjectStrand <- as.vector(subjectStrand)
            subjectRanges <- unname(ranges(subject))

            matchMatrix <-
              do.call(rbind,
                      lapply(commonSeqnames, function(seqnm)
                      {
                          qIdxs <- querySplitRanges[[seqnm]]
                          sIdxs <- subjectSplitRanges[[seqnm]]
                          overlaps <-
                            findOverlaps(seqselect(queryRanges, qIdxs),
                                         seqselect(subjectRanges, sIdxs),
                                         maxgap = maxgap, type = type,
                                         select = "all")
                          qHits <- queryHits(overlaps)
                          sHits <- subjectHits(overlaps)
                          matches <-
                            cbind(query = as.integer(qIdxs)[qHits],
                                  subject = as.integer(sIdxs)[sHits])
                          matches[which(seqselect(queryStrand, qIdxs)[qHits] *
                                        seqselect(subjectStrand, sIdxs)[sHits] != -1L), ,
                                  drop=FALSE]
                      }))
            matchMatrix <-
              matchMatrix[IRanges:::orderTwoIntegers(matchMatrix[ , 1L, drop=TRUE],
                                                     matchMatrix[ , 2L, drop=TRUE]), ,
                          drop=FALSE]
        }
        if (select == "all") {
            new("RangesMatching", matchMatrix = matchMatrix, DIM = DIM)
        } else {
            .matchMatrixToVector(matchMatrix, length(query))
        }
    }
)

setMethod("findOverlaps", c("GRangesList", "GRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        if (!IRanges:::isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)

        ans <-
          callGeneric(unlist(query, use.names=FALSE), subject,
                      maxgap = maxgap, type = type, select = "all")
        matchMatrix <- ans@matchMatrix
        matchMatrix[, 1L] <-
          togroup(query@partitioning)[matchMatrix[, 1L, drop=TRUE]]
        matchMatrix <- .cleanMatchMatrix(matchMatrix)
        if (select == "all") {
            DIM <- c(length(query), length(subject))
            initialize(ans, matchMatrix = matchMatrix, DIM = DIM)
        } else {
            .matchMatrixToVector(matchMatrix, length(query))
        }
    }
)

setMethod("findOverlaps", c("GRanges", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
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
        ans <-
          callGeneric(query, unlistSubject, maxgap = maxgap,
                      type = type, select = "all")
        matchMatrix <- ans@matchMatrix
        matchMatrix[, 2L] <- subjectGroups[matchMatrix[, 2L, drop=TRUE]]
        matchMatrix <- .cleanMatchMatrix(matchMatrix)
        if (select == "all") {
            DIM <- c(length(query), length(subject))
            initialize(ans, matchMatrix = matchMatrix, DIM = DIM)
        } else {
            .matchMatrixToVector(matchMatrix, length(query))
        }
    }
)

setMethod("findOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
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
        ans <-
          callGeneric(unlist(query, use.names=FALSE), unlistSubject,
                      maxgap = maxgap, type = type, select = "all")
        matchMatrix <- ans@matchMatrix
        matchMatrix[, 1L] <-
          togroup(query@partitioning)[matchMatrix[, 1L, drop=TRUE]]
        matchMatrix[, 2L] <- subjectGroups[matchMatrix[, 2L, drop=TRUE]]
        matchMatrix <- .cleanMatchMatrix(matchMatrix)
        if (select == "all") {
            DIM <- c(length(query), length(subject))
            initialize(ans, matchMatrix = matchMatrix, DIM = DIM)
        } else {
            .matchMatrixToVector(matchMatrix, length(query))
        }
    }
)


### =========================================================================
### countOverlaps methods
### -------------------------------------------------------------------------

setMethod("countOverlaps", c("GRanges", "GRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        type <- match.arg(type)
        tabulate(queryHits(findOverlaps(query, subject, maxgap = maxgap,
                                        type = type)), length(query))
    }
)

setMethod("countOverlaps", c("GRangesList", "GRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        type <- match.arg(type)
        tabulate(queryHits(findOverlaps(query, subject, maxgap = maxgap,
                                        type = type)), length(query))
    }
)

setMethod("countOverlaps", c("GRanges", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        type <- match.arg(type)
        tabulate(queryHits(findOverlaps(query, subject, maxgap = maxgap,
                                        type = type)), length(query))
    }
)

setMethod("countOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        type <- match.arg(type)
        tabulate(queryHits(findOverlaps(query, subject, maxgap = maxgap,
                                        type = type)), length(query))
    }
)


### =========================================================================
### subsetByOverlaps methods
### -------------------------------------------------------------------------

setMethod("subsetByOverlaps", c("GRanges", "GRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        type <- match.arg(type)
        query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                                  minoverlap = minoverlap, type = type,
                                  select = "first"))]
    }
)

setMethod("subsetByOverlaps", c("GRangesList", "GRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        type <- match.arg(type)
        query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                                  minoverlap = minoverlap, type = type,
                                  select = "first"))]
    }
)

setMethod("subsetByOverlaps", c("GRanges", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        type <- match.arg(type)
        query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                                  minoverlap = minoverlap, type = type,
                                  select = "first"))]
    }
)

setMethod("subsetByOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        type <- match.arg(type)
        query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                                  minoverlap = minoverlap, type = type,
                                  select = "first"))]
    }
)


### =========================================================================
### match methods
### -------------------------------------------------------------------------

setMethod("match", c("GRanges", "GRanges"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL)
    {
        if (length(nomatch) != 1)
            stop("'nomatch' must be of length 1") 
        ans <- findOverlaps(x, table, select = "first")
        if (!is.na(nomatch) && IRanges:::anyMissing(ans))
            ans[is.na(ans)] <- nomatch
        ans
    }
)

setMethod("match", c("GRanges", "GRangesList"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL)
    {
        if (length(nomatch) != 1)
            stop("'nomatch' must be of length 1") 
        ans <- findOverlaps(x, table, select = "first")
        if (!is.na(nomatch) && IRanges:::anyMissing(ans))
            ans[is.na(ans)] <- nomatch
        ans
    }
)

setMethod("match", c("GRangesList", "GRanges"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL)
    {
        if (length(nomatch) != 1)
            stop("'nomatch' must be of length 1") 
        ans <- findOverlaps(x, table, select = "first")
        if (!is.na(nomatch) && IRanges:::anyMissing(ans))
            ans[is.na(ans)] <- nomatch
        ans
    }
)

setMethod("match", c("GRangesList", "GRangesList"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL)
    {
        if (length(nomatch) != 1)
            stop("'nomatch' must be of length 1") 
        ans <- findOverlaps(x, table, select = "first")
        if (!is.na(nomatch) && IRanges:::anyMissing(ans))
            ans[is.na(ans)] <- nomatch
        ans
    }
)


### =========================================================================
### %in% methods
### -------------------------------------------------------------------------

setMethod("%in%", c("GRanges", "GRanges"),
    function(x, table)
        !is.na(match(x, table))
)

setMethod("%in%", c("GRanges", "GRangesList"),
    function(x, table)
        !is.na(match(x, table))
)

setMethod("%in%", c("GRangesList", "GRanges"),
    function(x, table)
        !is.na(match(x, table))
)

setMethod("%in%", c("GRangesList", "GRangesList"),
    function(x, table)
        !is.na(match(x, table))
)
