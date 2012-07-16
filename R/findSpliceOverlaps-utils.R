### =========================================================================
### findSpliceOverlaps utilities
### -------------------------------------------------------------------------

.rangeForSorted <- function(x) {
  part <- PartitioningByWidth(x)
  xflat <- unlist(x, use.names=FALSE)
  IRanges(start(xflat)[start(part)], end(xflat)[end(part)])
}

.compatibleTranscription <- function(query, splice, subject, intron, hits,
                                     clip=0L)
{
    ## FIXME : range() is slow
    qrange <- range(query)
    srange <- range(subject)
    splrange <- ranges(splice)
    qstrand <- as.character(strand(qrange))
    ## FIXME : should clip be a modifiable parameter?
    if (clip != 0L) {
        bounds <- .rangeForSorted(qrange) - clip
        qrange <- restrict(qrange, start(bounds), end(bounds),
                           keep.all.ranges = TRUE)
    }
    ## novelBounds
    qstrand <- as.character(qstrand)
    qrange <- unlist(qrange, use.names=FALSE)
    srange <- unlist(srange, use.names=FALSE)
    ## bounds violation for each match 
    lviol <- start(qrange) < start(srange)
    rviol <- end(qrange) > end(srange)
    TSSviol <- as.logical(ifelse(qstrand != "-", lviol, rviol))
    TSEviol <- as.logical(ifelse(qstrand != "+", lviol, rviol))
    ## bounds violation across all subjects hit (grouped)
    lmatch <- start(qrange) == start(srange)
    rmatch <- end(qrange) == end(srange)
    TSSmatch <- as.logical(ifelse(qstrand != "-", lmatch, rmatch))
    TSEmatch <- as.logical(ifelse(qstrand != "+", lmatch, rmatch))
    qhits <- queryHits(hits)
    TSSgroup <- (rowsum(as.integer(TSSmatch), qhits)[,1] > 0L)[factor(qhits)]
    TSEgroup <- (rowsum(as.integer(TSEmatch), qhits)[,1] > 0L)[factor(qhits)]

    TSS <- TSSviol & !TSSgroup
    TSE <- TSEviol & !TSEgroup
    boundsv <- TSSviol | TSEviol

    ## novelSplicing
    idx <- elementLengths(splrange) != 0L
    splicev <- rep.int(FALSE, length(splrange))
    if (any(idx)) {
        splrng <- splrange[idx]
        intrng <- ranges(intron)[idx]
        diff1 <- elementLengths(setdiff(splrng, intrng))
        diff2 <- elementLengths(setdiff(intrng, splrng))
        splicev[idx] <- (diff1 > 0L) | (diff2 > 0L)
    }

    DataFrame(compatible=!boundsv & !splicev, 
              novelTSS=TSS,
              novelTSE=TSE,
              novelSplicing=splicev)
}

.allMatch <- function(x, idx)
{
    (rowsum(as.integer(x), idx)[,1] == table(idx))[idx]
}

.oneMatch <- function(x, idx)
{
    xcnt <- rowsum(as.integer(x), idx)[,1]
    oneMatch <- rep((xcnt == 1L), table(idx))
    unname(x & oneMatch)
}

.novelExon <- function(splice, intron, nc)
{
    if (sum(elementLengths(splice)) == 0L)
        return(rep.int(FALSE, length(splice)))

    ## subset on elements with splices
    ans <- rep.int(FALSE, length(splice))
    idx <- elementLengths(splice) > 0 
    splice <- splice[idx]
    internal <- .gaps(splice)
    if (sum(elementLengths(internal)) == 0L)
        return(ans)

    hits <- findOverlaps(unlist(internal, use.names=FALSE), intron[idx], 
                         type="within", ignore.strand=TRUE)
    if (length(hits) > 0L) {
        ans0 <- rep.int(FALSE, length(splice))
        ne <- table(togroup(splice)[queryHits(hits)]) == 1L
        ans0[unique(togroup(splice)[queryHits(hits)])] <- ne
        ans[idx] <- ans0
        ans
    } else {
        ans
    }
}

.novelSpliceEvent <- function(splice, intron)
{
    if (sum(elementLengths(splice)) == 0L |
        sum(elementLengths(intron)) == 0L)
        DataFrame(Site=rep.int(FALSE, length(splice)),
                  Junction=rep.int(FALSE, length(splice)))

    ## subset on elements with splices
    site <- junction <- rep.int(FALSE, length(splice))
    idx <- elementLengths(splice) > 0 
    splice <- splice[idx]
    intron <- intron[idx]

    iflat <- unlist(intron, use.names=FALSE)
    sflat <- unlist(splice, use.names=FALSE)
    site[idx] <- .spliceEvent("site", iflat, sflat, splice)
    junction[idx] <- .spliceEvent("junction", iflat, sflat, splice)
    DataFrame(Site=site, Junction=junction)
}

.spliceEvent <- function(type, iflat, sflat, splice)
{
    if (type == "site") {
        combiner <- c
        elt <- rep(togroup(splice), each=2)
    } else if (type == "junction") {
        combiner <- paste
        elt <- togroup(splice)
    }
    ikeys <- paste(seqnames(iflat), combiner(start(iflat),
                   end(iflat)), strand(iflat), sep = ":")
    skeys <- paste(seqnames(sflat), combiner(start(sflat),
                   end(sflat)), strand(sflat), sep = ":")
    novel <- !(skeys %in% ikeys)
    rowsum(as.integer(novel), elt) > 0L
}

.novelRetention <- function(query, intron, ns)
{
    ## FIXME :  novelSplicing does not catch reads 
    ##          with no gaps

    ans <- rep.int(FALSE, length(query))
    if (sum(ns) == 0L)
        return(ans)

    ## subset on elements with novelSplicing
    qns <- query[ns]
    ins <- intron[ns] 
    ## FIXME : memory problem with paired-end
    hit <- findOverlaps(unlist(ins, use.names=FALSE), qns, 
                         type="within", ignore.strand=TRUE)
    if (length(hits) > 0L) {
        ans0 <- rep.int(FALSE, length(qns))
        nr <- table(togroup(qns)[queryHits(hit)]) == 1L
        ans0[unique(togroup(qns)[queryHits(hit)])] <- nr
        ans[ns] <- ans0
        ans
    } else {
        ans
    }
}

.result <- function(hits, nc=NULL, compatible=NULL, unique=NULL, coding=NULL, 
                    strandSpecific=NULL, novelTSS=NULL, novelTSE=NULL, 
                    novelSite=NULL, novelJunction=NULL, novelExon=NULL, 
                    novelRetention=NULL)
{
    nms <- c("compatible", "unique", "coding", "strandSpecific", 
             "novelTSS", "novelTSE", "novelSite", "novelJunction", 
             "novelExon", "novelRetention")
    ## full result
    if (!is.null(nc)) {
        values(hits) <- DataFrame(compatible, unique, coding, strandSpecific,
                                  novelTSS, novelTSE, novelSite, novelJunction, 
                                  novelExon, novelRetention) 
        hits
    ## no overlaps 
    } else if (is.null(compatible)) {
        mat <- matrix(logical(0), length(hits), length(nms))
        values(hits) <- DataFrame(mat)
        names(values(hits)) <- nms
        hits
    ## no compatible overlaps 
    } else {
        mat <- matrix(FALSE, length(hits), length(nms)) 
        values(hits) <- DataFrame(cbind(compatible, unique, coding, 
                                  strandSpecific, mat))
        names(values(hits)) <- nms
        hits
    }
}

## Until we have the formal 'gaps' method for GRangeList
isNumericOrNAs <- IRanges:::isNumericOrNAs
.gaps <- function(x, start=NA, end=NA)
{
    if (!isNumericOrNAs(start))
        stop("'start' must be an integer vector or NA")
    if (!is.integer(start))
        start <- as.integer(start)
    if (!isNumericOrNAs(end))
        stop("'end' must be an integer vector or NA")
    if (!is.integer(end))
        end <- as.integer(end)

    ## seqname and strand consistent in list elements
    if (all(elementLengths(runValue(seqnames(x))) == 1L) &&
        all(elementLengths(runValue(strand(x))) == 1L)) {
        flat <- unlist(x, use.names=FALSE)
        gaps <- gaps(ranges(x), start, end)

        idx <- elementLengths(gaps) != 0
        ## FIXME : can't handle lists with empty elements 
        ##         'start' and 'end' not quite right here
        firstseg <- start(PartitioningByWidth(x))
        seqnms <- rep(seqnames(flat)[firstseg], elementLengths(gaps))
        strand <- rep(strand(flat)[firstseg], elementLengths(gaps))
        gr <- relist(GRanges(seqnms, unlist(gaps, use.names=FALSE), strand), gaps)
        gr
    } else {
       psetdiff(range(x), x)
    }
}
