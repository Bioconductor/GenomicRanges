### =========================================================================
### Utility functions for performing basic ranges operations on GenomicRanges
### objects
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Intra-interval endomorphisms.
###

setMethod("shift", "GenomicRanges",
    function(x, shift=0L, use.names=TRUE)
    {
        ranges <- shift(ranges(x), shift, use.names=use.names)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            end(x, check=FALSE) <- end(ranges)
            start(x) <- pmin.int(start(ranges), end(x))
        } else {
            x <- clone(x, ranges=ranges)
        }
        x
    }
)

setMethod("narrow", "GenomicRanges",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        rng <- narrow(ranges(x), start=start, end=end, width=width,
                      use.names=TRUE)
        ranges(x) <- rng
        x
    }
)

setMethod("flank", "GenomicRanges",
    function(x, width, start=TRUE, both=FALSE, use.names=TRUE,
             ignore.strand=FALSE)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand)
            start <- rep.int(TRUE, length(x))
        else 
            start <- as.vector(start == (strand(x) != "-"))
        ranges <-
            flank(ranges(x), width=width, start=start, both=both,
                  use.names=use.names)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            start(x) <- start(ranges)
            end(x) <- end(ranges)
        } else {
            x <- clone(x, ranges=ranges)
        }
        x
    }
)

setMethod("resize", "GenomicRanges",
    function(x, width, fix="start", use.names=TRUE, ignore.strand=FALSE)
    {
        if (!missing(fix) &&
            (length(fix) > length(x) || length(x) %% length(fix) > 0L))
            stop("'x' is not a multiple of 'fix' length")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand) {
            fix <- Rle(rep.int(fix, length(x)))
        } else {
            revFix <- c(start="end", end="start", center="center")       
            fix <- ifelse(strand(x) == "-", revFix[fix], fix)
        }
        ranges <-
            resize(ranges(x), width=width, fix=fix, use.names=use.names)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            start(x) <- start(ranges)
            end(x) <- end(ranges)
        } else {
            x <- clone(x, ranges=ranges)
        }
        x
    }
)

## zooming (symmetrically scales the width)
setMethod("Ops", c("GenomicRanges", "numeric"),
    function(e1, e2)
    {
        if (IRanges:::anyMissing(e2))
            stop("NA not allowed as zoom factor")
        e2 <- recycleNumericArg(e2, "e2", length(e1))
        if (.Generic == "*") {
            e2 <- ifelse(e2 < 0, abs(1/e2), e2)
            resize(e1, width(e1) / e2, fix="center")
        } else {
            if (.Generic == "-") {
                e2 <- -e2
                .Generic <- "+"
            }
            if (.Generic == "+") {
                if (any(-e2*2 > width(e1)))
                    stop("adjustment would result in ranges ",
                         "with negative widths")
                resize(e1, width(e1) + e2*2, fix="center")
            }
        }
    }
)

.checkParms <- function(x, parm)
{
    if (!all(is.na(parm))) {
        if (!all(names(parm) %in% levels(seqnames(x))))
            stop("start should be a named numeric vector ",
                 "corresponding to seqnames")
    }
    temp <- structure(rep(NA_integer_, length(levels(seqnames(x)))), 
                      names=levels(seqnames(x)))
    temp[names(parm)] <- parm
    temp
}

.restrictRngs <- function(x, start, end, keep.all.ranges, use.names)
{
    tmp <- names(x)
    names(x) <- seq_len(length(x))
    rng <- ranges(x)
    res <- restrict(ranges(x), start, end, keep.all.ranges, use.names=TRUE)
    x <- x[as.numeric(names(res))]
    ranges(x) <- res
    if (!use.names)
        names(x) <- NULL
    else 
        names(x) <- tmp[as.numeric(names(res))]
    x
}

setMethod("restrict", "GenomicRanges",
    function(x, start=NA, end=NA, keep.all.ranges=FALSE, use.names=TRUE)
    {
        if (is.null(names(start)) && is.null(names(end)))
            return(.restrictRngs(x, start, end,keep.all.ranges, use.names))
        elementMetadata(x) <- cbind(elementMetadata(x),
                                    DataFrame(posIndx=seq_len(length(x))))
        splt <- split(x, seqnames(x))
        start <- .checkParms(x, start)
        end <- .checkParms(x, end) 
        res <- lapply(names(splt),
                      function(i) {
                          .restrictRngs(splt[[i]], start=start[i], end=end[i],
                                        keep.all.ranges, use.names)
                      })
        names(res) <- names(splt)
        ord <- unlist(GRangesList(res), use.names=FALSE)
        df <- data.frame(orig=elementMetadata(ord)$posIndx,
                         final=seq_len(length(ord)))
        indx <- with(df, order(orig, final))
        ord <- ord[indx, ]
        elementMetadata(ord) <- subset(elementMetadata(ord), select=-c(posIndx))
        ord
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Inter-interval endomorphisms.
###

.interIntervalGenomicRanges <- function(x, FUN, ignore.strand=FALSE, ...)
{
    x <- clone(x)
    elementMetadata(x) <- NULL
    if (ignore.strand)
       f <- paste(seqnames(x), Rle(factor("*"), length(x)), sep="\r")
    else
       f <- paste(seqnames(x), strand(x), sep="\r")
    xIRangesList <- split(unname(ranges(x)), f)
    ansIRangesList <- FUN(xIRangesList, ...)
    k <- elementLengths(ansIRangesList)
    splitListNames <- strsplit(names(ansIRangesList), split="\r", fixed=TRUE)
    listNameMatrix <- matrix(as.character(unlist(splitListNames)), nrow=2L)
    ansSeqnames <- Rle(factor(listNameMatrix[1L, ], levels=seqlevels(x)), k)
    ansStrand <- Rle(strand(listNameMatrix[2L, ]), k)
    update(x, seqnames=ansSeqnames,
           ranges=unlist(ansIRangesList, use.names=FALSE),
           strand=ansStrand,
           elementMetadata=new("DataFrame", nrows=length(ansSeqnames)))
}

setMethod("gaps", "GenomicRanges",
    function(x, start=1L, end=seqlengths(x))
    {
        seqlevels <- seqlevels(x)
        if (!is.null(names(start)))
            start <- start[seqlevels]
        if (!is.null(names(end)))
            end <- end[seqlevels]
        start <- IRanges:::recycleVector(start, length(seqlevels))
        start <- rep(start, each=3L)
        end <- IRanges:::recycleVector(end, length(seqlevels))
        end <- rep(end, each=3L)
        .interIntervalGenomicRanges(x, gaps, start=start, end=end)
    }
)

setMethod("range", "GenomicRanges",
    function(x, ..., ignore.strand=FALSE, na.rm=FALSE)
        .interIntervalGenomicRanges(unname(c(x, ...)), range, 
                                    ignore.strand=ignore.strand)
)

setMethod("reduce", "GenomicRanges",
    function(x, drop.empty.ranges=FALSE, min.gapwidth=1L,
             with.inframe.attrib=FALSE, ignore.strand=FALSE)
    {
        if (!identical(with.inframe.attrib, FALSE))
            stop("'with.inframe.attrib' argument not supported ",
                 "when reducing a GenomicRanges object")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand)
            strand(x) <- "*"
        .interIntervalGenomicRanges(x, reduce,
                                    drop.empty.ranges=drop.empty.ranges,
                                    min.gapwidth=min.gapwidth)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "isDisjoint", "disjoin", and "disjointBins" methods.
###

applyOnRangesBySpace <- function(x, FUN, ..., ignore.strand = FALSE) {
  if (ignore.strand)
    f <- seqnames(x)
  else
    f <- paste(seqnames(x), strand(x), sep = "\r")
  xIRangesList <- split(unname(ranges(x)), f)
  ans <- FUN(xIRangesList, ...)
  if (is(ans, "List")) # values per range, otherwise assumed to be per-space
    ans <- unsplit(ans, f)
  ans
}

setMethod("isDisjoint", "GenomicRanges",
    function(x, ignore.strand=FALSE)
    {
        all(applyOnRangesBySpace(x, isDisjoint, ignore.strand = ignore.strand))
    }
)

setMethod("disjoin", "GenomicRanges",
    function(x, ignore.strand=FALSE)
        .interIntervalGenomicRanges(x, disjoin, ignore.strand=ignore.strand)
)

setMethod("disjointBins", "GenomicRanges",
          function(x, ignore.strand = FALSE) {
            applyOnRangesBySpace(x, disjointBins,
                                 ignore.strand = ignore.strand)
          })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "precede", "follow", "nearest", "distance", and "distanceToNearest"
### methods.
###

.orderNumeric <- function(x)            # unstable order
    sort.list(x, na.last=NA, method="quick")

.GenomicRanges_findNearest0 <-
    function (query, subject, sentinel, leftOf=TRUE)
    ## query 'leftOf' subject
{
    sentinelidx <- length(subject) + seq_along(sentinel)
    subject <- c(subject, sentinel)

    ord <- .orderNumeric(subject)
    subject <- subject[ord]

    rle <- Rle(subject)
    subject <- runValue(rle)

    i <- findInterval(query - !leftOf, subject) + leftOf
    i[subject[i] %in% sentinel] <- NA_integer_
 
    IRanges:::.vectorToHits(i, rle, ord)
}

.findNearest_distance <-
    function(hit, query, subject, leftOf=TRUE)
    ## query 'leftOf' of subject
{
    if (leftOf)
        start(subject)[subjectHits(hit)] - end(query)[queryHits(hit)]
    else
        start(query)[queryHits(hit)] - end(subject)[subjectHits(hit)]
}

.Hits <-
    function(queryHits, subjectHits,
             queryLength=as.integer(max(c(0, queryHits))),
             subjectLength=as.integer(max(c(0, subjectHits))))
{
    o <- IRanges:::orderIntegerPairs(queryHits, subjectHits)
    new("Hits", queryHits=queryHits[o], subjectHits=subjectHits[o],
        queryLength=queryLength, subjectLength=subjectLength)
}

.findPrecedeFollow_pmin <-
    function(hit0, dist0, hit1, dist1)
    ## hit0, hit1 are (possibly multiple) hits on same query; repeated
    ## hits are equidistant within hit0 or hit1, but not necessarily
    ## between
{
    stopifnot(queryLength(hit0) == queryLength(hit1))
    stopifnot(subjectLength(hit0) == subjectLength(hit1))
    stopifnot(length(hit0) == length(dist0) ||
        length(hit1) == length(dist1))

    n <- queryLength(hit0)
    d0 <- d1 <- integer()
    d0[n] <- d1[n] <- NA_integer_

    i0 <- queryHits(hit0)
    i1 <- queryHits(hit1)

    d0[i0] <- dist0
    d1[i1] <- dist1

    dMin <- pmin.int(d0, d1, na.rm=TRUE)
    i0 <- dist0 == dMin[i0]
    i1 <- dist1 == dMin[i1]

    .Hits(c(queryHits(hit0)[i0], queryHits(hit1)[i1]),
        c(subjectHits(hit0)[i0], subjectHits(hit1)[i1]),
        queryLength(hit0), subjectLength(hit0))
}

.GenomicRanges_findPrecedeFollow <-
    function(query, subject, select, ignore.strand, 
             where=c("precede", "follow"))
{
    leftOf <- "precede" == match.arg(where)
    if (ignore.strand)
        strand(x) <- strand(subject) <- "+"

    if (leftOf) {
        plusfun <- function(xstart, xend, ystart, yend, sentinel)
            .GenomicRanges_findNearest0(xend, ystart, sentinel, leftOf)
        minusfun <- function(xstart, xend, ystart, yend, sentinel)
            .GenomicRanges_findNearest0(xstart, yend, sentinel, !leftOf)
    } else {
        plusfun <- function(xstart, xend, ystart, yend, sentinel)
            .GenomicRanges_findNearest0(xstart, yend, sentinel, leftOf)
        minusfun <- function(xstart, xend, ystart, yend, sentinel)
            .GenomicRanges_findNearest0(xend, ystart, sentinel, !leftOf)
    }

    ## sentinels marking seqlevels ends
    maxend <- max(max(end(query)), max(end(subject))) + 1
    lvls <- union(seqlevels(query), seqlevels(subject))
    offset <- setNames((seq_along(lvls) - 1) * maxend, lvls)
    stopifnot(typeof(offset) == "double")      # avoid integer overflow
    sentinel <- c(0, seq_along(lvls) * maxend)

    ## offset for sentinels
    queryOff <- unname(offset[as.character(seqnames(query))])
    queryStart <- start(query) + queryOff
    queryEnd <- end(query) + queryOff
    qid <- seq_along(query)

    subjectOff <- unname(offset[as.character(seqnames(subject))])
    subjectStart <- start(subject) + subjectOff
    subjectEnd <- end(subject) + subjectOff
    spid <- which(strand(subject) != "-")
    smid <- which(strand(subject) != "+")

    ## '+' query
    idx <- which(strand(query) == "+")
    phit <- plusfun(queryStart[idx], queryEnd[idx],
        subjectStart[spid], subjectEnd[spid], sentinel)
    phit <- .Hits(qid[idx][queryHits(phit)], spid[subjectHits(phit)])

    ## '-' query
    idx <- which(strand(query) == "-")
    mhit <- minusfun(queryStart[idx], queryEnd[idx],
        subjectStart[smid], subjectEnd[smid], sentinel)
    mhit <- .Hits(qid[idx][queryHits(mhit)], smid[subjectHits(mhit)])

    ## '*' query
    idx <- which(strand(query) == "*")
    bhit <- local({
        qid <- qid[idx]
        phit <- plusfun(queryStart[idx], queryEnd[idx],
            subjectStart[spid], subjectEnd[spid], sentinel)
        phit <- .Hits(qid[queryHits(phit)], spid[subjectHits(phit)],
            length(query), length(subject))
        mhit <- minusfun(queryStart[idx], queryEnd[idx],
            subjectStart[smid], subjectEnd[smid], sentinel)
        mhit <- .Hits(qid[queryHits(mhit)], smid[subjectHits(mhit)],
            length(query), length(subject))
        pdist <- .findNearest_distance(phit, query, subject, leftOf)
        mdist <- .findNearest_distance(mhit, query, subject, !leftOf)
        .findPrecedeFollow_pmin(phit, pdist, mhit, mdist)
    })

    ## clean up
    qryHits <- c(queryHits(phit), queryHits(mhit), queryHits(bhit))
    subjHits <- c(subjectHits(phit), subjectHits(mhit), subjectHits(bhit))
    if ("arbitrary" == select) {
        hits <- integer()
        hits[length(query)] <- NA_integer_
        idx <- !duplicated(qryHits) ## ties
        hits[qryHits[idx]] <- subjHits[idx]
    } else {
        hits <- .Hits(qryHits, subjHits, length(query), length(subject))
    }
    hits
}

setMethod("precede", c("GenomicRanges", "GenomicRanges"),
    function(x, subject, select = c("arbitrary", "all"), 
             ignore.strand=FALSE, ...)
    {
        select <- match.arg(select)
        .GenomicRanges_findPrecedeFollow(x, subject, select, ignore.strand, 
            "precede", ...) 
    }
)

setMethod("precede", c("GenomicRanges", "missing"),
    function(x, subject, select = c("arbitrary", "all"), 
             ignore.strand=FALSE, ...)
    {
        callGeneric(x, subject=x, select=select, ignore.strand=ignore.strand, 
            ...) 
    }
)

setMethod("follow", c("GenomicRanges", "GenomicRanges"),
    function(x, subject, select = c("arbitrary", "all"), 
             ignore.strand=FALSE, ...)
    {
        select <- match.arg(select)
        .GenomicRanges_findPrecedeFollow(x, subject, select, ignore.strand, 
            "follow", ...) 
    }
)

setMethod("follow", c("GenomicRanges", "missing"),
    function(x, subject, select=c("arbitrary", "all"), ignore.strand=FALSE, ...)
    {
        callGeneric(x, subject=x, select=select, ignore.strand=ignore.strand, 
            ...) 
    }
)

.nearest <- function(x, subject, ignore.strand, ...)
{
    p <- precede(x, subject, ignore.strand, select="arbitrary")
    f <- follow(x, subject, ignore.strand, select="arbitrary")
    midx <- !is.na(p) & !is.na(f)
    pdist <- distance(x[midx], subject[p[midx]])
    fdist <- distance(x[midx], subject[f[midx]])

    ## choose nearest or not missing 
    ans <- rep.int(NA_integer_, length(x))
    ans[midx] <- ifelse(pdist > fdist, p[midx], f[midx])
    ans[!midx] <- ifelse(is.na(p[!midx]), f[!midx], p[!midx])

    ## if tie, choose lowest index 
    didx <- pdist == fdist
    if (any(na.omit(didx)))
        ans[midx[didx]] <- pmin(p[midx[didx]], f[midx[didx]])

    ans
}

setMethod("nearest", c("GenomicRanges", "GenomicRanges"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        .nearest(x, subject, ignore.strand, ...) 
    }
)

setMethod("nearest", c("GenomicRanges", "missing"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        callGeneric(x, subject=x, ignore.strand=ignore.strand, ...)
    }
)

setMethod("distance", c("GenomicRanges", "GenomicRanges"),
    function(x, y, ignore.strand=FALSE, ...)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (length(x) != length(y))
            stop("'x' and 'y' must have the same length")
        d <- distance(ranges(x), ranges(y))
        mismtch <- as.character(seqnames(x)) != as.character(seqnames(y))
        if (any(mismtch))
            d[mismtch] <- NA
        if (!ignore.strand) {
            idx <- as.numeric(strand(x)) + as.numeric(strand(y))
            if (any(idx == 3))
                d[idx == 3] <- NA
        }
        d
    }
)

setMethod("distanceToNearest", c("GenomicRanges", "GenomicRanges"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        x_nearest <- nearest(x, subject, ignore.strand=ignore.strand)
        if (identical(integer(0), x_nearest)) {
            DataFrame(queryHits=integer(), subjectHits=integer(),
                      distance=integer())
        } else {
            queryHits=seq(length(x))[!is.na(x_nearest)]
            subjectHits=na.omit(x_nearest)
            distance=distance(x[!is.na(x_nearest)], subject[na.omit(x_nearest)],
                              ignore.strand=ignore.strand)
            DataFrame(queryHits, subjectHits, distance)
        }
    }
)

setMethod("distanceToNearest", c("GenomicRanges", "missing"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        callGeneric(x, subject=x, ignore.strand=ignore.strand, ...)
    }
)

