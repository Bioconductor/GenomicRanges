###

test_mergeNamedAtomicVectors <- function()
{
    ## The function is currently not exported.
    mergeNamedAtomicVectors <- GenomicRanges:::mergeNamedAtomicVectors

    x <- c(a=1, b=2, c=3, d=NA)
    checkIdentical(mergeNamedAtomicVectors(x, x), x)

    y1 <- x[FALSE]
    checkIdentical(mergeNamedAtomicVectors(x, y1), x)
    checkIdentical(mergeNamedAtomicVectors(y1, x), x)

    y2 <- c(c=NA_real_, d=NA_real_, b=NA_real_)
    checkIdentical(mergeNamedAtomicVectors(x, y2), x)
    checkIdentical(mergeNamedAtomicVectors(y2, x)[names(x)], x)
    checkIdentical(mergeNamedAtomicVectors(y2, x)[names(y2)], x[names(y2)])

    y3 <- c(c=3, e=NA, d=4, b=2)
    got <- mergeNamedAtomicVectors(x, y3)
    want <- c(a=1, b=2, c=3, d=4, e=NA)
    checkIdentical(got, want)
    got <- mergeNamedAtomicVectors(y3, x)
    want <- c(c=3, e=NA, d=4, b=2, a=1)
    checkIdentical(got, want)

    y4 <- c(c=0, e=NA, d=4, b=0)
    checkException(mergeNamedAtomicVectors(x, y4), silent=TRUE)
    checkException(mergeNamedAtomicVectors(y4, x), silent=TRUE)

    x2 <- c(a=1, b=NA, c=3, d=NA)
    y5 <- c(c=NA, e=5, d=4, b=2)
    got <- mergeNamedAtomicVectors(x2, y5)
    want <- c(a=1, b=2, c=3, d=4, e=5)
    checkIdentical(got, want)
}

test_keepSeqlevels <- function()
{
    gr <- GRanges(c("chr1", "chr1", "chr2", "chr3"), IRanges(1:4, width=3))
    ## unnamed
    checkIdentical("chr1", seqlevels(keepSeqlevels(gr, "chr1")))
    got <- seqlevels(keepSeqlevels(gr, c("chr1", "chr3")))
    checkIdentical(c("chr1", "chr3"), got)
    ## named
    checkIdentical("chr1", seqlevels(keepSeqlevels(gr, c(foo="chr1"))))
    ## bogus
    got <- seqlevels(suppressWarnings(keepSeqlevels(gr, "chrX")))
    checkIdentical(character(0), got)
}

test_dropSeqlevels <- function()
{
    gr <- GRanges(c("chr1", "chr1", "chr2", "chr3"), IRanges(1:4, width=3))
    ## unnamed
    checkIdentical(c("chr2", "chr3"), seqlevels(dropSeqlevels(gr, "chr1")))
    got <- seqlevels(dropSeqlevels(gr, c("chr1", "chr3")))
    checkIdentical("chr2", got)
    ## named
    got <- seqlevels(dropSeqlevels(gr, c(foo="chr1")))
    checkIdentical(c("chr2", "chr3"), got)
    ## bogus
    got <- seqlevels(suppressWarnings(dropSeqlevels(gr, "chrX")))
    checkIdentical(seqlevels(gr), got)

    grl <- split(gr, as.character(seqnames(gr)))
    got <- dropSeqlevels(grl, c("chr1", "chr2"))
    checkIdentical("chr3", seqlevels(got))
    checkIdentical(1L, length(got))
}

test_renameSeqlevels <- function()
{
    gr <- GRanges(c("chr1", "chr1", "chr2", "chr3"), IRanges(1:4, width=3))
    checkException(renameSeqlevels(gr, "CHR1"), silent=TRUE)
    got <- seqlevels(renameSeqlevels(gr, c("chr2", "CHR1", "chr3")))
    checkIdentical(c("chr2", "CHR1", "chr3"), got)
    got <- seqlevels(suppressWarnings(renameSeqlevels(gr, c(foo="chr2"))))
    checkIdentical(seqlevels(gr), got)
    got <- seqlevels(renameSeqlevels(gr, c(chr2="CHR2")))
    checkIdentical(c("chr1", "CHR2", "chr3"), got)
}

