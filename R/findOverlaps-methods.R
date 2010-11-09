## No longer needed because we are no longer second guessing the users intent.
## .similarSeqnameConvention <- function(seqs1, seqs2) {
##     funList <-
##       list(isRoman = function(x) grepl("[ivIV]$", x),
##            isArabic = function(x) grepl("[0-9]$", x),
##            haschr = function(x) grepl("^chr", x),
##            hasChr = function(x) grepl("^Chr", x),
##            hasCHR = function(x) grepl("^CHR", x),
##            haschrom = function(x) grepl("^chrom", x),
##            hasChrom = function(x) grepl("^Chrom", x),
##            hasCHROM = function(x) grepl("^CHROM", x))
##     all(sapply(funList, function(f) any(f(seqs1)) == any(f(seqs2))))
## }

.testSeqEquiv <- function(seqs1, seqs2) {
  length(setdiff(seqs1, seqs2))!=0 &&
  length(setdiff(seqs2, seqs1))!=0
}                                         


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
        stop("overlap type \"", type, "\" is unsupported ",
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
    overlaps01 <-findOverlaps(query0, subject1,
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
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        #        if (!identical(minoverlap, 1L))
        #    warning("'minoverlap' argument is ignored")
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
            ## TODO: Use merge() when it becomes available.
            #seqinfo <- merge(seqinfo(query), seqinfo(subject))
            seqinfo <- seqinfo(query)
            querySeqnames <- seqnames(query)
            querySplitRanges <- splitRanges(querySeqnames)
            uniqueQuerySeqnames <- names(querySplitRanges)

            subjectSeqnames <- seqnames(subject)
            subjectSplitRanges <- splitRanges(subjectSeqnames)
            uniqueSubjectSeqnames <- names(subjectSplitRanges)

            ## uniqueQuerySeqnames & uniqueSubjectSeqnames are my sets
            ## It's ok if both sets are equivalent.
            ## Also ok if one set is contained within the other.
            ## But not ok if there are elements from each set that are not
            ## in the other set.
            if(.testSeqEquiv(uniqueQuerySeqnames, uniqueSubjectSeqnames)){ 
                ## A different warning if there are no matches at all
                if(is.na(table(uniqueQuerySeqnames
                               %in% uniqueSubjectSeqnames)["TRUE"])){
                  warning("no seqnames from 'query' and 'subject' were identical")
                }else ## Versus having some things match...
                {
                  warning("some seqnames from 'query' and 'subject' differ")
                }
            }
            
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

setMethod("findOverlaps", c("GRangesList", "GenomicRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        if (!IRanges:::isSingleNumber(maxgap) || maxgap < 0)
            stop("'maxgap' must be a non-negative integer")
        type <- match.arg(type)
        select <- match.arg(select)
        unlistQuery <- unlist(query, use.names = FALSE)
        queryGroups <- togroup(query@partitioning)
        ans <-
          callGeneric(unlistQuery, subject,
                      maxgap = maxgap, type = type, select = "all")
        
        matchMatrix <- ans@matchMatrix
        if(minoverlap > 1L && nrow(ans@matchMatrix) > 0) {  
            matchMatrix <- ans@matchMatrix[FALSE,]
            intrsct <- pintersect(ranges(unlistQuery)[queryHits(ans)],
                    ranges(subject[subjectHits(ans)]))

            df <- data.frame(q = queryHits(ans),
                         subject = subjectHits(ans),
                         query =queryGroups[queryHits(ans)],
                         w=width(intrsct))
        
            mat <-  with(df, aggregate( w, list(query = query, subject = subject), sum))
            indx <- mat$x >= minoverlap
            if(any(indx)) 
                matchMatrix  <- 
                    as.matrix(mat[indx, c("query", "subject"), drop = FALSE], 
                            rownames.force = FALSE)
        } else {
            matchMatrix[, 1L] <-
                togroup(query@partitioning)[matchMatrix[, 1L, drop=TRUE]]
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

setMethod("findOverlaps", c("GenomicRanges", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
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
        ans <-
          callGeneric(query, unlistSubject, maxgap = maxgap,
                      type = type, select = "all")
        matchMatrix <- ans@matchMatrix
        if(minoverlap > 1L && nrow(ans@matchMatrix) > 0) {    
            matchMatrix <- ans@matchMatrix[FALSE,]
            intrsct <- pintersect(ranges(query)[queryHits(ans)],
                        ranges(unlistSubject[subjectHits(ans)]))
            df <- data.frame(query = queryHits(ans),
                         s = subjectHits(ans),
                         subject =subjectGroups[subjectHits(ans)],
                         w=width(intrsct))
            mat <-  with(df, aggregate( w, list(query = query, subject = subject), sum))
            indx <- mat$x >= minoverlap
            if(any(indx)) 
                matchMatrix  <- 
                    as.matrix(mat[indx, c("query", "subject"), drop = FALSE], 
                            rownames.force = FALSE)
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

setMethod("findOverlaps", c("GRangesList", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
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
            unlistSubject <-  unlistSubject[keep]
            subjectGroups <- subjectGroups[keep]
        } else if (type == "end") {
            keep <- end(subject@partitioning)[elementLengths(subject) > 0L]
            unlistSubject <-  unlistSubject[keep]
            subjectGroups <- subjectGroups[keep]
        }
        
        ans <-
          callGeneric(unlistQuery, unlistSubject,
                      maxgap = maxgap, type = type, select = "all")
        
        matchMatrix <- ans@matchMatrix
        if(minoverlap  > 1L  && nrow(ans@matchMatrix) > 0) {
            matchMatrix <- ans@matchMatrix[FALSE,]
            intrsct <- pintersect(ranges(unlistQuery)[queryHits(ans)],
                    ranges(unlistSubject[subjectHits(ans)]))

            df <- data.frame(q = queryHits(ans),
                         s = subjectHits(ans),
                         query = queryGroups[queryHits(ans)],
                         subject =subjectGroups[subjectHits(ans)],
                         w=width(intrsct))
            mat <-  with(df, aggregate( w, list(query = query, subject = subject), sum))
            indx <- mat$x >= minoverlap
            if(any(indx)) 
                matchMatrix  <- 
                    as.matrix(mat[indx, c("query", "subject"), drop = FALSE], 
                            rownames.force = FALSE)
        } else {
            matchMatrix[, 1L] <-
                togroup(query@partitioning)[matchMatrix[, 1L, drop=TRUE]]
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

setMethod("findOverlaps", c("RangesList", "GenomicRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        findOverlaps(as(query, "GRanges"), subject = subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select))
    }
)

setMethod("findOverlaps", c("RangesList", "GRangesList"),
          function(query, subject, maxgap = 0L, minoverlap = 1L,
                   type = c("any", "start", "end"),
                   select = c("all", "first"))
          {
            findOverlaps(as(query, "GRanges"), subject = subject,
                         maxgap = maxgap, minoverlap = minoverlap,
                         type = match.arg(type), select = match.arg(select))
          })

setMethod("findOverlaps", c("GenomicRanges", "RangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        findOverlaps(query, as(subject, "GRanges"),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select))
    }
)

setMethod("findOverlaps", c("GRangesList", "RangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        findOverlaps(query, as(subject, "GRanges"),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select))
    }
)

setMethod("findOverlaps", c("RangedData", "GenomicRanges"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        ## Calls "findOverlaps" method for c("RangesList", "GenomicRanges")
        ## defined above.
        findOverlaps(ranges(query), subject = subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select))
    }
)

setMethod("findOverlaps", c("RangedData", "GRangesList"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        ## Calls "findOverlaps" method for c("RangesList", "GRangesList")
        ## defined above.
        findOverlaps(ranges(query), subject = subject,
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select))
    }
)

setMethod("findOverlaps", c("GenomicRanges", "RangedData"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        ## Calls "findOverlaps" method for c("GenomicRanges", "RangesList")
        ## defined above.
        findOverlaps(query, ranges(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select))
    }
)

setMethod("findOverlaps", c("GRangesList", "RangedData"),
    function(query, subject, maxgap = 0L, minoverlap = 1L,
             type = c("any", "start", "end"),
             select = c("all", "first"))
    {
        ## Calls "findOverlaps" method for c("GRangesList", "RangesList")
        ## defined above.
        findOverlaps(query, ranges(subject),
                     maxgap = maxgap, minoverlap = minoverlap,
                     type = match.arg(type), select = match.arg(select))
    }
)


### =========================================================================
### findOverlaps-based methods
### -------------------------------------------------------------------------

.countOverlaps.default <- function(query, subject,
                                   maxgap = 0L, minoverlap = 1L,
                                   type = c("any", "start", "end"))
{
    if (!identical(minoverlap, 1L))
        warning("'minoverlap' argument is ignored")
    type <- match.arg(type)
    tabulate(queryHits(findOverlaps(query, subject, maxgap = maxgap,
                                    type = type)), length(query))
}

.subsetByOverlaps.default <- function(query, subject,
                                      maxgap = 0L, minoverlap = 1L,
                                      type = c("any", "start", "end"))
{
    if (!identical(minoverlap, 1L))
        warning("'minoverlap' argument is ignored")
    type <- match.arg(type)
    query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                              minoverlap = minoverlap, type = type,
                              select = "first"))]
}

.subsetByOverlaps2 <- function(query, subject, maxgap = 0L, minoverlap = 1L,
                               type = c("any", "start", "end"))
{
    if (!identical(minoverlap, 1L))
        warning("'minoverlap' argument is ignored")
    type <- match.arg(type)
    i <- !is.na(findOverlaps(query, subject, maxgap = maxgap,
                             minoverlap = minoverlap, type = type,
                             select = "first"))
    query[seqsplit(i, space(query), drop=FALSE)]
}

.subsetByOverlaps3 <- function(query, subject, maxgap = 0L, minoverlap = 1L,
                               type = c("any", "start", "end"))
{
  if (!identical(minoverlap, 1L))
    warning("'minoverlap' argument is ignored")
  type <- match.arg(type)
  query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                            minoverlap = minoverlap, type = type,
                            select = "first")),]
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
    c("GRangesList", "RangedData")
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

