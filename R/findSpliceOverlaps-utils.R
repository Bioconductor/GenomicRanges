### =========================================================================
### findSpliceOverlaps utilities
### -------------------------------------------------------------------------

.readRanges <- function(bam, param, singleEnd)
{
    if (!"XS" %in% bamTag(param))
        bamTag(param) <- c(bamTag(param), "XS")
    if (singleEnd)
        reads <- readGappedAlignments(path(bam), param=param)
    else
        reads <- readGappedAlignmentPairs(path(bam), param=param)

    metadata(reads)$bamfile <- bam
    ## adjust strand based on 'XS'
    if (!is.null(xs <- values(reads)$XS))
        strand(reads) <- ifelse(!is.na(xs) & xs != "?", xs, "*")
    reads
}

.rangeForSorted <- function(x) {
  part <- PartitioningByWidth(x)
  x_flat <- unlist(x, use.names=FALSE)
  IRanges(start(x_flat)[start(part)], end(x_flat)[end(part)])
}

.compatibleTranscription <- function(query, subject, hits, clip=0L)
{
    qrange <- ranges(query)
    srange <- ranges(subject)
    splice <- ranges(.gaps(query))
    intron <- ranges(.gaps(subject))
    if (clip != 0L) {
        bounds <- .rangeForSorted(qrange) - clip
        qrange <- restrict(qrange, start(bounds), end(bounds),
                           keep.all.ranges = TRUE)
    }
    boundsv <- elementLengths(setdiff(qrange, srange)) != 0L
    idx <- elementLengths(splice) != 0L
    splicev <- rep.int(FALSE, length(splice))
    if (any(idx))
        splicev[idx] <- elementLengths(setdiff(splice[idx], intron[idx])) > 0L

    DataFrame(compatible=!boundsv & !splicev, 
              novelBounds=boundsv,
              novelSplicing=splicev)
}

.allMatch <- function(x, idx)
{
    (rowsum(as.integer(x), idx)[,1] == table(idx))[idx]
}

.oneMatch <- function(x, idx)
{
    xcnt <- rowsum(as.integer(x), idx)[,1]
    oneMatch <- (xcnt == 1L)[idx]
    unname(x & oneMatch)
}

.novelBounds <- function(query, subject, hits)
## Assuming no geneID available in annotation
## Subjects are grouped according to the query they hit
{
    ## FIXME : novelBounds for performance?
    q <- unlist(range(query), use.names=FALSE)
    s <- unlist(range(subject), use.names=FALSE)

    ## bounds violation for each match 
    lviol <- start(q) < start(s)
    rviol <- end(q) > end(s)
    TSSviol <- as.logical(ifelse(strand(q) == "+", lviol, rviol))
    TSEviol <- as.logical(ifelse(strand(q) == "-", lviol, rviol))
    ## bounds violation across all subjects hit (grouped)
    lmatch <- start(q) == start(s)
    rmatch <- end(q) == end(s)
    TSSmatch <- as.logical(ifelse(strand(q) == "+", lmatch, rmatch))
    TSEmatch <- as.logical(ifelse(strand(q) == "-", lmatch, rmatch))
    qhits <- queryHits(hits)
    TSSgroup <- (rowsum(as.integer(TSSmatch), qhits)[,1] > 0L)[qhits] 
    TSEgroup <- (rowsum(as.integer(TSEmatch), qhits)[,1] > 0L)[qhits]

    DataFrame(TSS=TSSviol & !TSSgroup, TSE=TSEviol & !TSEgroup)
}

.novelExon <- function(splice, intron, nc)
{
    ## FIXME : splice[nc] or novelSplicing for performance?
    if (sum(elementLengths(splice)) == 0L)
        return(rep.int(FALSE, length(splice)))
    internal <- .gaps(splice)
    ans <- rep.int(FALSE, length(splice))
    if (sum(elementLengths(internal)) == 0L)
        return(ans)

    hits <- findOverlaps(unlist(internal, use.names=FALSE), intron, 
                         type="within", ignore.strand=TRUE)
    if (length(hits) > 0L) {
        ne <- table(togroup(splice)[queryHits(hits)]) == 1L
        ans[unique(togroup(splice)[queryHits(hits)])] <- ne
        ans
    } else {
        ans
    }
}

.novelRetention <- function(query, intron, nc)
{
    ## FIXME : query[nc] or novelSplicing for performance?
    hits <- findOverlaps(unlist(intron, use.names=FALSE), query, 
                         type="within", ignore.strand=TRUE)
    ans <- rep.int(FALSE, length(query))
    if (length(hits) > 0L) {
        nr <- table(togroup(query)[queryHits(hits)]) == 1L
        ans[unique(togroup(query)[queryHits(hits)])] <- nr
        ans
    } else {
        ans
    }
}

.novelSpliceEvent <- function(splice, intron)
{
    ## FIXME : novelSplicing for performance?
    if (sum(elementLengths(splice)) == 0L |
        sum(elementLengths(intron)) == 0L) {
        DataFrame(Site=rep.int(FALSE, length(splice)),
                  Junction=rep.int(FALSE, length(splice)))
    } else {
        iflat <- unlist(intron, use.names=FALSE)
        sflat <- unlist(splice, use.names=FALSE)
        DataFrame(Site=.spliceEvent("site", iflat, sflat, splice),
                  Junction=.spliceEvent("junction", iflat, sflat, splice))
    }
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

.result <- function(hits, nc=NULL, compatible=NULL, unique=NULL, coding=NULL, 
                    strandSpecific=NULL, novelSplicing=NULL, 
                    novelTSS=NULL, novelTSE=NULL, novelSite=NULL, 
                    novelJunction=NULL, novelExon=NULL, novelRetention=NULL)
{
    nms <- c("compatible", "unique", "coding", "strandSpecific", 
             "novelSplicing", "novelTSS", "novelTSE", "novelSite", 
             "novelJunction", "novelExon", "novelRetention")
    ## full result
    if (!is.null(nc)) {
        values(hits) <- DataFrame(compatible, unique, coding, strandSpecific,
                                  novelSplicing, novelTSS, novelTSE, novelSite,
                                  novelJunction, novelExon, novelRetention) 
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
        firstseg <- start(PartitioningByWidth(x))
        seqnms <- rep(seqnames(flat)[firstseg], elementLengths(gaps))
        strand <- rep(strand(flat)[firstseg], elementLengths(gaps))
        relist(GRanges(seqnms, unlist(gaps, use.names=FALSE), strand), gaps)
    } else {
        psetdiff(range(x), x)
    }
}
