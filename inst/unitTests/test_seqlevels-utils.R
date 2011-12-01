gr1 <- GRanges(seqnames = c("chr1", "chr2"),
          ranges = IRanges(start=1, width = 3))
gr2 <- GRanges(seqnames = c("chr1", "chr1", "chr2", "chr3", "chr3"),
          ranges = IRanges(start=c(10, 20, 1, 1, 10), width=3))

test_keepSeqlevels <- function()
{
    grl <- GRangesList(GRanges("A", IRanges(1, 1)))
    metadata(grl) <- list(x=1)
    checkIdentical(metadata(grl), metadata(keepSeqlevels(grl, "A"))) 

    grl <- GRangesList(gr1, gr2)
    checkIdentical(3L, length(unlist(keepSeqlevels(grl, "chr1"))))

    gr <- GRanges(seqnames = c("chr1", "chr2"),
                  ranges = IRanges(start=1, width=3))
    checkIdentical(seqlengths(gr)[1], seqlengths(keepSeqlevels(gr, "chr1")))

    galn <- GappedAlignments(names = c("A","B"), rname = Rle(c("chr1", "chr2")),
                             pos = as.integer(c(10, 100)), cigar = c("50M", "50M"),
                             strand=strand(c("*", "*")))
    metadata(galn) <- list(x=1)
    checkIdentical(metadata(galn), metadata(keepSeqlevels(galn, "chr1"))) 
}

test_renameSeqlevels <- function()
{
    grl <- GRangesList(GRanges("A", IRanges(1, 1)))
    metadata(grl) <- list(x=1)
    checkIdentical(metadata(grl), metadata(renameSeqlevels(grl, c(A="a"))))

    grl <- GRangesList(gr1, gr2)
    checkIdentical(c("CHR1", "CHR2", "CHR3"), seqlevels(renameSeqlevels(grl, 
        c(chr1="CHR1", chr2="CHR2", chr3="CHR3"))))

    checkIdentical(c("chr1", "chr1"), 
        seqlevels(suppressWarnings(renameSeqlevels(gr1, c(chr2="chr1")))))
}

