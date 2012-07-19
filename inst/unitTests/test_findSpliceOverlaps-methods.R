.extract <- function(x, col) as.logical(values(x)[[col]])

test_findSpliceOverlaps_novelTSS <- function()
{
    ## strand
    genes <- GRangesList(
        GRanges("chr1", IRanges(5, 15), "+"),
        GRanges("chr1", IRanges(5, 15), "-"),
        GRanges("chr1", IRanges(5, 15), "*"))
    reads <- GRangesList(GRanges("chr1", IRanges(3, 13), "+"))
    res <- findSpliceOverlaps(reads, genes[1])
    checkIdentical(TRUE, .extract(res, "novelTSS"))
    res <- suppressWarnings(findSpliceOverlaps(reads, genes[2]))
    checkIdentical(logical(0), .extract(res, "novelTSS"))
    res <- findSpliceOverlaps(reads, genes[3])
    checkIdentical(TRUE, .extract(res, "novelTSS"))

    ## multiple matches 
    genes <- GRangesList(
        GRanges("chr1", IRanges(15, 20), "+"),
        GRanges("chr1", IRanges(10, 20), "+"),
        GRanges("chr1", IRanges(5, 20), "+"))
    reads <- GRangesList(GRanges("chr1", IRanges(5, 20), "+"))
    res <- findSpliceOverlaps(reads, genes[1])
    checkIdentical(TRUE, .extract(res, "novelTSS"))
    res <- findSpliceOverlaps(reads, genes[1:2])
    checkIdentical(c(TRUE, TRUE), .extract(res, "novelTSS"))
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(c(FALSE, FALSE, FALSE), .extract(res, "novelTSS"))

    ## gaps 
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 15), c(10, 20)), "+"))
    reads <- GRangesList(
        GRanges("chr1", IRanges(12, 23), "+"),
        GRanges("chr1", IRanges(3, 18), "+"))
    res <- findSpliceOverlaps(reads[1], genes)
    checkIdentical(FALSE, .extract(res, "novelTSS"))
    res <- findSpliceOverlaps(reads[2], genes)
    checkIdentical(TRUE, .extract(res, "novelTSS"))
}

test_findSpliceOverlaps_novelTSE <- function()
{
    ## strand
    genes <- GRangesList(
        GRanges("chr1", IRanges(5, 15), "+"),
        GRanges("chr1", IRanges(5, 15), "-"),
        GRanges("chr1", IRanges(5, 15), "*"))
    reads <- GRangesList(GRanges("chr1", IRanges(12, 18), "+"))
    res <- findSpliceOverlaps(reads, genes[1])
    checkIdentical(TRUE, .extract(res, "novelTSE"))
    res <- findSpliceOverlaps(reads, genes[2])
    checkIdentical(logical(0), .extract(res, "novelTSE"))
    res <- findSpliceOverlaps(reads, genes[3])
    checkIdentical(TRUE, .extract(res, "novelTSE"))

    ## multiple matches 
    genes <- GRangesList(
        GRanges("chr1", IRanges(5, 15), "+"),
        GRanges("chr1", IRanges(5, 20), "+"),
        GRanges("chr1", IRanges(5, 25), "+"))
    reads <- GRangesList(GRanges("chr1", IRanges(5, 25), "+"))
    res <- findSpliceOverlaps(reads, genes[1])
    checkIdentical(TRUE, .extract(res, "novelTSE"))
    res <- findSpliceOverlaps(reads, genes[1:2])
    checkIdentical(c(TRUE, TRUE), .extract(res, "novelTSE"))
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(c(FALSE, FALSE, FALSE), .extract(res, "novelTSE"))

    ## gaps 
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 15), c(10, 20)), "+"))
    reads <- GRangesList(
        GRanges("chr1", IRanges(2, 12), "+"),
        GRanges("chr1", IRanges(18, 23), "+"))
    res <- findSpliceOverlaps(reads[1], genes)
    checkIdentical(FALSE, .extract(res, "novelTSE"))
    res <- findSpliceOverlaps(reads[2], genes)
    checkIdentical(TRUE, .extract(res, "novelTSE"))
}

test_findSpliceOverlaps_novelExon <- function()
{
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"))
    ## 'within' intron boundaries
    reads <- GRangesList(
        GRanges("chr1", IRanges(c(7, 12, 20), c(10, 18, 23)), "+"))
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(TRUE, .extract(res, "novelExon"))

    ## FIXME :currently TRUE
    ## Do we want a novel exon to be completely w/in?
    reads <- GRangesList(
        GRanges("chr1", IRanges(c(5, 9, 20), c(7, 12, 23)), "+"),
        GRanges("chr1", IRanges(c(7, 15, 23), c(10, 21, 25)), "+"))
    res <- findSpliceOverlaps(reads[1], genes)
    checkIdentical(TRUE, .extract(res, "novelExon"))
    res <- findSpliceOverlaps(reads[2], genes)
    checkIdentical(TRUE, .extract(res, "novelExon"))

    ## region not 'intronic' in all transcripts
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"),
        GRanges("chr1", IRanges(c(5, 20), c(15, 25)), "+"))
    reads <- GRangesList(
        GRanges("chr1", IRanges(c(7, 12, 20), c(10, 18, 23)), "+"))
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(c(FALSE, FALSE), .extract(res, "novelExon"))
}

test_findSpliceOverlaps_novelRetention <- function()
{
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"))
    ## 'within' intron boundaries
    reads <- GRangesList(
        GRanges("chr1", IRanges(c(7, 12, 20), c(10, 18, 23)), "+"))
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(TRUE, .extract(res, "novelRetention"))

    ## exact match to intron boundaries
    ## (no overlap between 'reads' and 'genes')
    reads <- GRangesList(
        GRanges("chr1", IRanges(11, 19), "+"))
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(logical(0), .extract(res, "compatible"))

    ## overlap and span intron boundaries
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 20, 30), c(10, 25, 35)), "+"))
    reads <- GRangesList(
        GRanges("chr1", IRanges(5, 12), "+"), 
        GRanges("chr1", IRanges(18, 23), "+"),
        GRanges("chr1", IRanges(c(4, 30), c(26, 36)), "+"))
    res <- findSpliceOverlaps(reads[1], genes) ## no read gaps
    checkIdentical(TRUE, .extract(res, "novelRetention"))
    res <- findSpliceOverlaps(reads[3], genes) ## read gaps
    checkIdentical(TRUE, .extract(res, "novelRetention"))

    ## FIXME : hits a portion of the intronic region
    ##         but is not completely 'within'
    ## region is not 'intronic' in all transcripts  
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"),
        GRanges("chr1", IRanges(c(5, 20), c(15, 25)), "+"))
    reads <- GRangesList(
        GRanges("chr1", IRanges(4, 26), "+"))
    res <- findSpliceOverlaps(reads[1], genes)
    checkIdentical(c(TRUE, TRUE), .extract(res, "novelRetention"))
}

test_findSpliceOverlaps_novelSite <- function()
{
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 15), c(10, 20)), "+"))
    ## single novel site, novel junction
    reads <- GRangesList(
        GRanges("chr1", IRanges(c(5, 15), c(7, 20)), "+"))
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(TRUE, .extract(res, "novelSite"))
    checkIdentical(TRUE, .extract(res, "novelJunction"))

    ## two novel sites, novel junction
    reads <- GRangesList(
        GRanges("chr1", IRanges(c(5, 17), c(7, 20)), "+"))
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(TRUE, .extract(res, "novelSite"))
    checkIdentical(TRUE, .extract(res, "novelJunction"))
}

test_findSpliceOverlaps_novelJunction <- function()
{
    ## novel junction, no novel sites
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"),
        GRanges("chr1", IRanges(c(5, 22), c(15, 25)), "+"))

    ## query = GRanges 
    reads <- GRangesList(
        GRanges("chr1", IRanges(c(5, 20), c(15, 25)), "+"))
    GRres <- findSpliceOverlaps(reads, genes)
    checkIdentical(c(TRUE, TRUE), .extract(GRres, "novelJunction"))
    checkIdentical(c(FALSE, FALSE), .extract(GRres, "novelSite"))

    ## query = GappedAlignments
    gal <- GappedAlignments("chr1", 5L, "11M4N6M", strand("+"))
    GALres <- findSpliceOverlaps(gal, genes)
    checkIdentical(c(TRUE, TRUE), .extract(GALres, "novelJunction"))
    checkIdentical(c(FALSE, FALSE), .extract(GALres, "novelSite"))

    ## query = GappedAlignmentPairs
    gal1 <- GappedAlignments("chr1", 5L, "11M4N6M", strand("+"))
    gal2 <- GappedAlignments("chr1", 50L, "6M", strand("-"))
    galp <- GappedAlignmentPairs(gal1, gal2, TRUE)
    GALPres <- findSpliceOverlaps(galp, genes)
    checkIdentical(c(TRUE, TRUE), .extract(GALPres, "novelJunction"))
    checkIdentical(c(FALSE, FALSE), .extract(GALPres, "novelSite"))
}

