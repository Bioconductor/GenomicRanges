### =========================================================================
### nearest (and related) methods
### -------------------------------------------------------------------------
###
### Dependencies :
###
###        distanceToNearest
###          |          |
###       nearest    distance
###       |     |
###  precede   follow

.orderNumeric <- function(x)            # unstable order
    sort.list(x, na.last=NA, method="quick")

.GenomicRanges_findNearest0 <-
    function (query, subject, sentinel, leftOf=TRUE)
    ## query 'leftOf' subject
{
    zeroInSubject <- subject %in% 0L
    if (any(zeroInSubject))
        subject <- unique(c(subject, sentinel))
    else
        subject <- c(subject, sentinel)
    ord <- .orderNumeric(subject)
    subject <- subject[ord]

    rle <- Rle(subject)
    subject <- runValue(rle)

    ## zero in query
    i <- findInterval(query - !leftOf, subject) + leftOf
    zeroInQuery <- i %in% 0L
    if (any(zeroInQuery)) {
        if (leftOf)
            i[zeroInQuery] <- 2L
        else
            i[zeroInQuery] <- 1L
    }
    ## zero in subject
    if (any(zeroInSubject))
        i[subject[i] %in% sentinel & subject[i] != 0L] <- NA_integer_
    else
        i[subject[i] %in% sentinel] <- NA_integer_

    IRanges:::vectorToHits(i, rle, ord)
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

.Hits <- function(queryHits, subjectHits,
                  queryLength=as.integer(max(c(0, queryHits))),
                  subjectLength=as.integer(max(c(0, subjectHits))))
{
    o <- orderIntegerPairs(queryHits, subjectHits)
    Hits(queryHits[o], subjectHits[o], queryLength, subjectLength,
         sort.by.query=TRUE)
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
    if (!length(query) || !length(subject))
        return(Hits(nLnode=length(query), nRnode=length(subject),
                    sort.by.query=TRUE))

    leftOf <- "precede" == match.arg(where)
    if (ignore.strand)
        strand(query) <- strand(subject) <- "+"

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
    endq <- start(query) + width(query) # end(query) incorrect for 0-width
    ends <- start(subject) + width(subject)
    maxend <- max(max(endq), max(ends)) + 1
    lvls <- union(seqlevels(query), seqlevels(subject))
    offset <- setNames((seq_along(lvls) - 1) * maxend, lvls)
    stopifnot(typeof(offset) == "double")      # avoid integer overflow
    sentinel <- c(0, seq_along(lvls) * maxend)

    ## offset for sentinels
    queryOff <- unname(offset[as.character(seqnames(query))])
    queryStart <- start(query) + queryOff
    queryEnd <- end(query) + queryOff   # true end + offset
    qid <- seq_along(query)

    subjectOff <- unname(offset[as.character(seqnames(subject))])
    subjectStart <- start(subject) + subjectOff
    subjectEnd <- end(subject) + subjectOff # true end + offset
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
    starphit <- plusfun(queryStart[idx], queryEnd[idx],
                    subjectStart[spid], subjectEnd[spid], sentinel)
    starphit <- .Hits(qid[idx][queryHits(starphit)],
                      spid[subjectHits(starphit)])
    ## '*' query and '*' subject pairs are treated as if on '+' strand;
    ## omit '*' subjects from this test
    smid <- which(strand(subject) == "-")
    starmhit <- minusfun(queryStart[idx], queryEnd[idx],
                     subjectStart[smid], subjectEnd[smid], sentinel)
    starmhit <- .Hits(qid[idx][queryHits(starmhit)], smid[subjectHits(starmhit)])
    starhit <- .Hits(c(queryHits(starphit), queryHits(starmhit)),
                     c(subjectHits(starphit), subjectHits(starmhit)),
                     length(query), length(subject))

    ## '*' strand query can return a value for both mhit and phit.
    ## Choose the closest range regardless of strand.
    both <- (idx %in% queryHits(starmhit)) & (idx %in% queryHits(starphit))

    if (any(both)) {
        x <- query[idx[both]]
        sidx <- which(queryHits(starhit) %in% idx[both])
        y <- subject[subjectHits(starhit)[sidx]]
        repeats <- tabulate(queryHits(starhit))[idx[both]]
        dist <- distance(rep(x, times=repeats), y)
        dist_il <- relist(dist, PartitioningByWidth(repeats))
        sidx_il <- relist(sidx, PartitioningByWidth(repeats))
        drops_il <- sidx_il[!dist_il == min(dist_il)]
        drop <- unlist(drops_il)
        if (length(drop))
            starhit <- starhit[-drop]
    }

    qryHits <- c(queryHits(phit), queryHits(mhit), queryHits(starhit))
    subjHits <- c(subjectHits(phit), subjectHits(mhit), subjectHits(starhit))
    hits <- .Hits(qryHits, subjHits, length(query), length(subject))
    ## Break ties
    if (select == "all") {
        hits
    } else if (select == "first") {
        first <- rep(NA_integer_, length(query))
        idx <- which(!duplicated(queryHits(hits)))
        first[queryHits(hits)[idx]] <- subjectHits(hits)[idx]
        first
    } else if ("last" == select) {
        last <- rep(NA_integer_, length(query))
        rev_query <- rev(queryHits(hits)) ## can't call rev() on Hits
        idx <- which(!duplicated(rev_query))
        last[rev_query[idx]] <- rev(subjectHits(hits))[idx]
        last
    } else
        stop("'select' must be one of c('first', 'last', 'all')")
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### precede() and follow()
###

setMethod("precede", c("GenomicRanges", "GenomicRanges"),
    function(x, subject,
             select=c("first", "all"),
             ignore.strand=FALSE)
    {
        select <- match.arg(select)
        .GenomicRanges_findPrecedeFollow(x, subject, select, ignore.strand,
                                         "precede")
    }
)

setMethod("precede", c("GenomicRanges", "missing"),
    function(x, subject,
             select=c("first", "all"),
             ignore.strand=FALSE)
    {
        select <- match.arg(select)
        .GenomicRanges_findPrecedeFollow(x, subject, select, ignore.strand,
                                         "precede")
    }
)

setMethod("follow", c("GenomicRanges", "GenomicRanges"),
    function(x, subject,
             select=c("last", "all"),
             ignore.strand=FALSE)
    {
        select <- match.arg(select)
        .GenomicRanges_findPrecedeFollow(x, subject, select, ignore.strand,
                                         "follow")
    }
)

setMethod("follow", c("GenomicRanges", "missing"),
    function(x, subject,
             select=c("last", "all"),
             ignore.strand=FALSE)
    {
        select <- match.arg(select)
        .GenomicRanges_findPrecedeFollow(x, subject, select, ignore.strand,
                                         "follow")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### nearest()
###

.filterHits <- function(hits, i, map)
{
    m <- as.matrix(hits[as(hits, "IRanges")[i]])
    m[, 1L] <- map[m[, 1L]]
    m
}

.nearest <- function(x, subject, select, ignore.strand, drop.self=FALSE)
{
    ## overlapping ranges
    if (drop.self) {
        ol <- findOverlaps(x, maxgap=0L, select=select,
                           ignore.strand=ignore.strand, drop.self=TRUE)
    } else {
        ol <- findOverlaps(x, subject, maxgap=0L, select=select,
                           ignore.strand=ignore.strand)
    }

    if (select == "all") {
        olv <- selectHits(ol, select="first")
    } else {
        olv <- ol
    }

    ## non-overlapping ranges
    if (length(x <- x[is.na(olv)]) != 0) {
        ## precede() and follow() do not support select="arbitrary"
        if (select == "arbitrary") {
            p <- precede(x, subject, select="first", ignore.strand)
            f <- follow(x, subject, select="last", ignore.strand)
        } else {
            p <- precede(x, subject, select, ignore.strand)
            f <- follow(x, subject, select, ignore.strand)
        }

        ## terminate if no results
        if (!length(p) && !length(f)) {
            if (is(olv, "Hits") && !length(olv) || all(is.na(olv))) {
                if (select == "all")
                    return(Hits(nLnode=length(x),
                                nRnode=length(subject),
                                sort.by.query=TRUE))
                else if (select == "arbitrary")
                    return (rep(NA, length(x)))
            }
        }

        if (select == "all") {
            p0 <- p
            p <- selectHits(p, select="first")
            f0 <- f
            f <- selectHits(f, select="last")
        }

        ## choose nearest or not missing
        pdist <- .nearestDistance(x, subject, p)
        fdist <- .nearestDistance(x, subject, f)
        pnearest <- ifelse(pdist == fdist, p < f, pdist < fdist)
        isNA <- is.na(pnearest)
        pnearest[isNA] <- is.na(f[isNA])

        if (select == "all") {
            map <- which(is.na(olv))
            pnearest[pdist == fdist] <- TRUE
            m <- rbind(as.matrix(ol), .filterHits(p0, pnearest, map),
                                      .filterHits(f0, !pnearest, map))
            m <- m[orderIntegerPairs(m[, 1L], m[, 2L]),, drop=FALSE]
            ol@from <- unname(m[, 1L])
            ol@to <- unname(m[, 2L])
        } else {
            olv[is.na(olv)] <- ifelse(pnearest, p, f)
            ol <- olv
        }
    }
    ol
}

.nearestDistance <- function(x, subject, index)
{
    if (length(index)) {
        maxStart <- pmax.int(start(x), start(subject)[index])
        minEnd <- pmin.int(end(x), end(subject)[index])
        pmax.int(maxStart - minEnd - 1L, 0L)
    } else NA
}

setMethod("nearest", c("GenomicRanges", "GenomicRanges"),
    function(x, subject,
             select=c("arbitrary", "all"),
             ignore.strand=FALSE)
    {
        select <- match.arg(select)
        .nearest(x, subject, select=select, ignore.strand=ignore.strand)
    }
)

setMethod("nearest", c("GenomicRanges", "missing"),
    function(x, subject,
             select=c("arbitrary", "all"),
             ignore.strand=FALSE)
    {
        select <- match.arg(select)
        .nearest(x, x, select=select, ignore.strand=ignore.strand,
                 drop.self=TRUE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### distance()
###

setMethod("distance", c("GenomicRanges", "GenomicRanges"),
    function(x, y, ignore.strand=FALSE, ...)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")

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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### distanceToNearest()
###

setMethod("distanceToNearest", c("GenomicRanges", "GenomicRanges"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        x_nearest <- nearest(x, subject, ignore.strand=ignore.strand, ...)
        .distanceToNearest(x_nearest, x, subject, ignore.strand=ignore.strand)
    }
)

setMethod("distanceToNearest", c("GenomicRanges", "missing"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        x_nearest <- nearest(x, ignore.strand=ignore.strand, ...)
        .distanceToNearest(x_nearest, x, x, ignore.strand=ignore.strand)
    }
)

.distanceToNearest <- function(x_nearest, x, subject, ignore.strand)
{
    ## 'x_nearest' is Hits when select = all
    if (is(x_nearest, "Hits")) {
        queryHits <- queryHits(x_nearest)
        subjectHits <- subjectHits(x_nearest)
    } else {
    ## 'x_nearest' is Integer vector when select = arbitrary
        queryHits <- seq_along(x)[!is.na(x_nearest)]
        subjectHits <- x_nearest[!is.na(x_nearest)]
    }

    if (!length(subjectHits) || all(is.na(subjectHits))) {
        Hits(nLnode=length(x),
             nRnode=length(subject),
             distance=integer(0),
             sort.by.query=TRUE)
    } else {
        distance <- distance(x[queryHits], subject[subjectHits],
                             ignore.strand=ignore.strand)
        Hits(queryHits, subjectHits, length(x), length(subject), distance,
             sort.by.query=TRUE)
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### precedes() and follows()
###

.normBounds <- function(x, si) {
    if (is(x, "GRangesList")) {
        x <- range(x)
        if (all(lengths(x) == 1L)) {
            x <- unlist(x)
        } else {
            stop("operation undefined when ranges cross seqnames and strands")
        }
    }
    seqlevels(x) <- seqlevels(si)
    x
}

precedes <- function(x, y) {
    si <- merge(seqinfo(x), seqinfo(y))
    x <- .normBounds(x, si)
    y <- .normBounds(y, si)
    seqnames(x) == seqnames(y) &
        ifelse(strand(y) == "-", start(x) > end(y), end(x) < start(y))
}

follows <- function(x, y) {
    si <- merge(seqinfo(x), seqinfo(y))
    x <- .normBounds(x, si)
    y <- .normBounds(y, si)
    seqnames(x) == seqnames(y) &
        ifelse(strand(y) == "-", end(x) < start(y), start(x) > end(y))
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Find 'k' nearest neighbors
###

.nearestKNeighbors_all_k <- function(x, k) {
    ## the 'k' nearest neighbors in 'x', when all equally-distant
    ## neighbors increment k by 1. Each element of 'x' is ordered. E.g.,
    ##
    ##   x = IntegerList(c(1, 1, 2), 1:2)
    ##   k = 1
    ##   return: c(2L, 1L)
    ##
    ans <- integer(length(x))
    non_zero_lengths <- lengths(x) > 0
    if (k == 0L || !any(non_zero_lengths))
        return(ans)

    x <- x[non_zero_lengths]
    idx <- pmin(k, lengths(x))
    value <- x[as.list(idx)]
    ans[non_zero_lengths] <- sum(x <= value)
    ans
}

.nearestKNeighbors <- function(x, subject, k, select, ignore.strand, drop.self=FALSE)
{
    seqlevels(subject) <- seqlevels(x)

    starts <- with(subject, GRanges(seqnames, IRanges(start, width=1L), strand))
    ends <- with(subject, GRanges(seqnames, IRanges(end, width=1L), strand))

    if (ignore.strand) {
        starts <- unstrand(starts)
        ends <- unstrand(ends)
    }

    start_ord <- order(starts)
    end_ord <- order(ends)

    starts <- starts[start_ord]
    ends <- ends[end_ord]

    ## NOTE: select="all" is needed in general for nearestKNeighbors here
    if (drop.self) {
        ol <- findOverlaps(x, maxgap=0L, select="all",
                           ignore.strand=ignore.strand, drop.self=TRUE)
    } else {
        ol <- findOverlaps(x, subject, maxgap=0L, select="all",
                           ignore.strand=ignore.strand)
    }
    ol_hits <- countLnodeHits(ol)
    ol <- as(ol, "List")
    overlaps <- IntegerList(lapply(ol_hits, rep, x=0L))

    if (length(x)) {
        phits <- precede(x, starts, ignore.strand=ignore.strand)
        fhits <- follow(x, ends, ignore.strand=ignore.strand)
    } else {
        phits <- fhits <- integer()
    }

    if (!ignore.strand) {
        exchange <- decode(strand(x) == "-")
        tmp <- phits[exchange]
        phits[exchange] <- fhits[exchange]
        fhits[exchange] <- tmp
    }

    findPart <- function(x, w) {
        S4Vectors:::findIntervalAndStartFromWidth(x, w)[["interval"]]
    }

    b <- width(disjoin(c(ranges(seqnames(starts)), ranges(strand(starts)))))
    width <- ifelse(select == "all", max(length(subject) - 1L, 1L), k)

    seqends <- end(strand(starts))[findPart(phits, b)]
    phits[is.na(phits)] <- 1L
    seqends[is.na(seqends)] <- 0L
    pwindows <- restrict(IRanges(phits, width = width), end=seqends)
    pwindows_kept <- seqends >= phits

    seqstarts <- start(strand(ends))[findPart(fhits, b)]
    seqstarts[is.na(seqstarts)] <- 1L
    fhits[is.na(fhits)] <- 0L
    fwindows <- restrict(IRanges(end=fhits, width = width), seqstarts)
    fwindows_kept <- seqstarts <= fhits

    pdist <- extractList(start(starts), pwindows) - end(x)
    pdist[!pwindows_kept] <- IntegerList(integer(0))
    fdist <- start(x) - extractList(end(ends), fwindows)
    fdist[!fwindows_kept] <- IntegerList(integer(0))

    dist <- pc(pdist, fdist, overlaps)
    dist <- abs(dist)

    pans <- extractList(start_ord, pwindows)
    pans[!pwindows_kept] <- IntegerList(integer(0))
    fans <- extractList(end_ord, fwindows)
    fans[!fwindows_kept] <- IntegerList(integer(0))

    ans <- pc(pans, fans, ol)

    o <- order(dist)
    if (select == "all")
        k <- .nearestKNeighbors_all_k(dist[o], k)
    ans <- ans[heads(o, k)]
    ans[lengths(ans) == 0L] <- NA
    ans
}

setGeneric("nearestKNeighbors", function(x, subject, ...) standardGeneric("nearestKNeighbors"))

setMethod("nearestKNeighbors", c("GenomicRanges", "missing"),
    function(x, subject, k = 1L, select = c("arbitrary", "all"),
             ignore.strand = FALSE)
{
    select <- match.arg(select)
    .nearestKNeighbors(x, x, k = k, select = select, ignore.strand = ignore.strand,
             drop.self = TRUE)
})

setMethod("nearestKNeighbors", c("GenomicRanges", "GenomicRanges"),
    function(x, subject, k = 1L, select = c("arbitrary", "all"),
             ignore.strand = FALSE)
{
    select <- match.arg(select)
    .nearestKNeighbors(x, subject, k = k, select = select, ignore.strand = ignore.strand)
})
