make_test_GRanges <- function() {
    new("GRanges",
        seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")),
        elementMetadata = DataFrame(score = 1:10, GC = seq(1, 0, length=10)))
}

test_GenomicRanges_inter_interval_ops <- function()
{
    ## gaps
    gr <- unname(make_test_GRanges())[ , character(0)]
    checkIdentical(gaps(gr, start = 1, end = 10),
                   GRanges(seqnames = Rle(c("chr1", "chr2", "chr3"), c(2, 3, 3)),
                           ranges = IRanges(start=1,
                                            end=c(5, 4, 1, 10, 3, 6, 8, 10)),
                           strand =
                           strand(c("+", "*", "+", "-", "*", "+", "-", "*"))))

    ## range
    gr <- unname(make_test_GRanges())[ , character(0)]
    checkIdentical(range(gr),
                   GRanges(seqnames = Rle(c("chr1", "chr2", "chr3"), c(3, 2, 2)),
                           ranges = IRanges(start=c(6, 1, 5, 2, 4, 7, 9), end=10),
                           strand = strand(c("+", "-", "*", "+", "*", "+", "-"))))
    checkTrue(validObject(range(gr, ignore.strand=TRUE)))

    ## reduce
    gr <- unname(make_test_GRanges())[ , character(0)]
    checkIdentical(reduce(gr),
                   GRanges(seqnames = Rle(c("chr1", "chr2", "chr3"), c(3, 2, 2)),
                           ranges = IRanges(start=c(6, 1, 5, 2, 4, 7, 9), end=10),
                           strand = strand(c("+", "-", "*", "+", "*", "+", "-"))))

    ## disjoin
    gr <- unname(make_test_GRanges())[ , character(0)]
    checkIdentical(disjoin(gr),
                   GRanges(seqnames = Rle(c("chr1", "chr2", "chr3"), c(3, 3, 4)),
                           ranges = IRanges(start=c(6, 1, 5, 2, 3, 4, 7, 8, 9, 10),
                                            end=c(10, 10, 10, 2, 10, 10, 7, 10, 9, 10)),
                           strand =
                           strand(c("+", "-", "*", "+", "+", "*", "+", "+", "-", "-"))))
}

