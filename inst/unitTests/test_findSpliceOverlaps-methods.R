findSpliceOverlaps <- GenomicRanges:::.findSpliceOverlaps 
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
        GRanges("chr1", IRanges(12, 18), "+"),
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
        GRanges("chr1", IRanges(18, 25), "+"))
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

    ## exact match to intron boundaries
    ## (FALSE b/c 'reads' has no gaps)
    reads <- GRangesList(
        GRanges("chr1", IRanges(c(7, 11, 20), c(10, 19, 23)), "+"))
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(FALSE, .extract(res, "novelExon"))

    ## overlap intron boundaries
    reads <- GRangesList(
        GRanges("chr1", IRanges(c(5, 9, 20), c(7, 12, 23)), "+"),
        GRanges("chr1", IRanges(c(7, 15, 23), c(10, 21, 25)), "+"))
    res <- findSpliceOverlaps(reads[1], genes)
    checkIdentical(FALSE, .extract(res, "novelExon"))
    res <- findSpliceOverlaps(reads[2], genes)
    checkIdentical(FALSE, .extract(res, "novelExon"))

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
    checkIdentical(FALSE, .extract(res, "novelRetention"))

    ## exact match to intron boundaries
    ## (no overlap between 'reads' and 'genes')
    reads <- GRangesList(
        GRanges("chr1", IRanges(11, 19), "+"))
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(logical(0), .extract(res, "novelRetention"))

    ## overlap and span intron boundaries
    reads <- GRangesList(
        GRanges("chr1", IRanges(5, 12), "+"),
        GRanges("chr1", IRanges(18, 23), "+"),
        GRanges("chr1", IRanges(4, 26), "+"))
    res <- findSpliceOverlaps(reads[1], genes)
    checkIdentical(FALSE, .extract(res, "novelRetention"))
    res <- findSpliceOverlaps(reads[2], genes)
    checkIdentical(FALSE, .extract(res, "novelRetention"))
    res <- findSpliceOverlaps(reads[3], genes)
    checkIdentical(TRUE, .extract(res, "novelRetention"))

    ## region is not 'intronic' in all transcripts  
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"),
        GRanges("chr1", IRanges(c(5, 20), c(15, 25)), "+"))
    reads <- GRangesList(
        GRanges("chr1", IRanges(4, 26), "+"))
    res <- findSpliceOverlaps(reads[1], genes)
    checkIdentical(c(FALSE, FALSE), .extract(res, "novelRetention"))
}

test_findSpliceOverlaps_novelSplicing <- function()
{
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 20, 35), c(10, 25, 40)), "+"))
    reads <- GRangesList(
        GRanges("chr1", IRanges(c(5, 24), c(22, 40)), "+"),
        GRanges("chr1", IRanges(c(5, 35), c(25, 40)), "+"))
    res <- findSpliceOverlaps(reads[1], genes)
    checkIdentical(TRUE, .extract(res, "novelSplicing"))
    res <- findSpliceOverlaps(reads[2], genes)
    checkIdentical(FALSE, .extract(res, "novelSplicing"))

    ## read with no gaps
    reads <- GRangesList(
        GRanges("chr1", IRanges(5, 40), "+"))
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(FALSE, .extract(res, "novelSplicing"))
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
    genes <- GRangesList(
        GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"),
        GRanges("chr1", IRanges(c(5, 22), c(15, 25)), "+"))
    reads <- GRangesList(
        GRanges("chr1", IRanges(c(5, 20), c(15, 25)), "+"))
    ## novel junction, no novel sites
    res <- findSpliceOverlaps(reads, genes)
    checkIdentical(c(TRUE, TRUE), .extract(res, "novelJunction"))
    checkIdentical(c(FALSE, FALSE), .extract(res, "novelSite"))
}

