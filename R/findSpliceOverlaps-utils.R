### =========================================================================
### findSpliceOverlaps utilities
### -------------------------------------------------------------------------

.compatibleTranscription <- function(query, subject, splice, clip=0L)
{
    qrng <- ranges(query)
    srng <- ranges(subject)
    sprng <- ranges(splice)
    ## FIXME : should clip be a modifiable parameter?
    ## Or we could consider taking it out entirely... the aligner should clip
    if (clip != 0L) {
        bounds <- .rangeForSorted(qrng) - clip
        qrange <- restrict(qrng, start(bounds), end(bounds),
                           keep.all.ranges = TRUE)
    }
    bnds <- elementLengths(setdiff(qrng, srng)) == 0L
    splc <- elementLengths(intersect(srng, sprng)) == 0L
    bnds & splc
}

.novelBounds <- function(query, subject, qhits)
{
    qrange <- range(query)
    qstrand <- as.character(strand(qrange))
    qrange <- unlist(qrange, use.names=FALSE)
    srange <- unlist(range(subject), use.names=FALSE)
    ## bounds violation for each match 
    lviol <- start(qrange) < start(srange)
    rviol <- end(qrange) > end(srange)
    TSSviol <- as.logical(ifelse(qstrand != "-", lviol, rviol))
    TSEviol <- as.logical(ifelse(qstrand != "+", lviol, rviol))
    ## bounds violation across all subjects hit (grouped)
    lviol <- start(qrange) == start(srange)
    rviol <- end(qrange) == end(srange)
    TSSmatch <- as.logical(ifelse(qstrand != "-", lviol, rviol))
    TSEmatch <- as.logical(ifelse(qstrand != "+", lviol, rviol))
    TSSgroup <- (rowsum(as.integer(TSSmatch), qhits)[,1] > 0L)[factor(qhits)]
    TSEgroup <- (rowsum(as.integer(TSEmatch), qhits)[,1] > 0L)[factor(qhits)]

    TSS <- TSSviol & !TSSgroup
    TSE <- TSEviol & !TSEgroup

    DataFrame(TSS=TSS, TSE=TSE)
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

.novelExon <- function(splice, intronRegion)
{
    if (sum(elementLengths(splice)) == 0L)
        return(logical(length(splice)))

    ## subset on elements with splices
    ans <- logical(length(splice))
    idx <- elementLengths(splice) > 0 
    splice <- splice[idx]
    internal <- unlist(.gaps(splice), use.names=FALSE)
    if (sum(length(internal)) == 0L)
        return(ans)

    ## FIXME : not competely "within"
    hits <- findOverlaps(internal, intronRegion, ignore.strand=TRUE)
    if (length(hits) > 0L) {
        ans0 <- logical(length(splice))
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

.novelRetention <- function(query, intronRegion)
{
    ans <- logical(length(query))
    if (length(query) == 0L)
        return(ans)

    hits <- findOverlaps(unlist(query, use.names=FALSE), intronRegion, 
                        ignore.strand=TRUE)
    if (length(hits) > 0L) {
        ans[unique(togroup(query)[queryHits(hits)])] <- TRUE
    } else {
        ans
    }
}

.intronicRegions <- function(tx, intron) {
  txflt <- unlist(tx, use.names = FALSE)
  intronflt <- unlist(intron, use.names = FALSE)
  regions <- setdiff(intronflt, txflt, ignore.strand = TRUE)
  #map <- findOverlaps(regions, intronflt)
  #mcols(regions)$tx_id <-
  #  seqsplit(names(tx)[togroup(introns)][subjectHits(intronic_to_tx)],
  #           queryHits(intronic_to_tx))
  regions
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
        mcols(hits) <- DataFrame(compatible, unique, coding, strandSpecific,
                                  novelTSS, novelTSE, novelSite, novelJunction,
                                  novelExon, novelRetention)
        hits
    ## no overlaps 
    } else if (is.null(compatible)) {
        mat <- matrix(logical(0), length(hits), length(nms))
        mcols(hits) <- DataFrame(mat)
        names(mcols(hits)) <- nms
        hits
    ## no compatible overlaps 
    } else {
        mat <- matrix(FALSE, length(hits), length(nms))
        mcols(hits) <- DataFrame(cbind(compatible, unique, coding,
                                  strandSpecific, mat))
        names(mcols(hits)) <- nms
        hits
    }
}

insertGaps <- function(reads) {
  query.break <- mcols(reads)$query.break
  if (is.null(query.break))
    stop("missing 'query.break' metadata variable: reads not paired?")
  reads_flat <- unlist(reads, use.names = FALSE)
  reads_part <- PartitioningByWidth(reads)
  left_end <- start(reads_part) + query.break - 1L
  right_start <- left_end + 1L
  start <- end(reads_flat)[left_end]
  end <- pmax(start(reads_flat)[right_start], start - 1L)
  if (any(seqnames(reads_flat)[left_end] != seqnames(reads_flat)[right_start]))
    stop("reads are on different chromosomes")
  GRanges(seqnames(reads_flat)[left_end],
          IRanges(start, end), strand(reads_flat)[left_end])
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
### FIXME: this makes this function more of an 'introns' than a .gaps.
### FIXME: this breaks when the GRangesList is not ordered by position
        if (!is.null(mcols(x)$query.break)) {
          insert_gaps <- split(ranges(insertGaps(x)), seq_len(length(x)))
          gaps <- setdiff(gaps, insert_gaps)
        }
        
        idx <- elementLengths(gaps) != 0
        ## FIXME : can't handle lists with empty elements 
        ##         'start' and 'end' not quite right here
        firstseg <- start(PartitioningByWidth(x))
        seqnms <- rep(seqnames(flat)[firstseg], elementLengths(gaps))
        strand <- rep(strand(flat)[firstseg], elementLengths(gaps))
        gr <- relist(GRanges(seqnms, unlist(gaps, use.names=FALSE), strand), gaps)
        gr
    } else {
### FIXME: does not handle query.break column yet
        psetdiff(range(x), x)
    }
    
}

.rangeForSorted <- function(x) {
  part <- PartitioningByWidth(x)
  xflat <- unlist(x, use.names=FALSE)
  IRanges(start(xflat)[start(part)], end(xflat)[end(part)])
}
