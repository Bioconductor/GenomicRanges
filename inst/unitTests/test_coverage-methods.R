make_test_GRanges <- function()
    GRanges(Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
            IRanges(1:10, end=10, names=head(letters, 10)),
            Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
            seqinfo=Seqinfo(paste0("chr", 1:3)),
            score=1:10, GC=seq(1, 0, length=10))

make_test_GRangesList <- function() {
    a <- make_test_GRanges()
    b <- GRanges(Rle(factor(c("chr2", "chr4", "chr5")), c(3, 6, 4)),
                 IRanges(1:13, end=13, names=tail(letters, 13)),
                 Rle(strand(c("-", "+", "-")), c(4, 5, 4)),
                 seqinfo=Seqinfo(paste0("chr", c(2, 4:5))),
                 score=1:13, GC=seq(0, 1, length=13))
    GRangesList(a=a, b=b)
}

test_coverage_GRanges <- function() {
    gr <- make_test_GRanges()

    target <- RleList(chr1=Rle(1:3, c(4, 1, 5)),
                      chr2=Rle(0:3, c(1, 1, 1, 7)),
                      chr3=Rle(0:4, c(6, 1, 1, 1, 1)),
                      compress=FALSE)
    checkIdentical(target, coverage(gr))

    width <- c(chr1=10, chr2=20, chr3=30)
    target <- RleList(chr1=Rle(1:3, c(4, 1, 5)),
                      chr2=Rle(c(0:3, 0L), c(1, 1, 1, 7, 10)),
                      chr3=Rle(c(0:4, 0L), c(6, 1, 1, 1, 1, 20)),
                      compress=FALSE)
    checkIdentical(target, coverage(gr, width=width))

    weight <- list(chr1=1L, chr2=10L, chr3=100L)
    target <- RleList(chr1=Rle(1:3, c(4, 1, 5)),
                      chr2=Rle(10L * 0:3, c(1, 1, 1, 7)),
                      chr3=Rle(100L * 0:4, c(6, 1, 1, 1, 1)),
                      compress=FALSE)
    checkIdentical(target, coverage(gr, weight=weight))

    shift <- list(chr1=0, chr2=1, chr3=2)
    target <- RleList(chr1=Rle(1:3, c(4, 1, 5)),
                      chr2=Rle(0:3, c(2, 1, 1, 7)),
                      chr3=Rle(0:4, c(8, 1, 1, 1, 1)),
                      compress=FALSE)
    checkIdentical(target, coverage(gr, shift=shift))

    ## with circular sequences
    gr <- GRanges(seqnames=c("A", "B"),
                  ranges=IRanges(start=5:6, width=7))
    gr@seqinfo <- Seqinfo(seqnames=c("A", "B"),
                          seqlengths=c(10, NA),
                          isCircular=c(TRUE, FALSE))
    target <- RleList(A=Rle(c(1L, 0L, 1L), c(1, 3, 6)),
                      B=Rle(c(0L, 1L), c(5, 7)),
                      compress=FALSE)
    checkIdentical(target, coverage(gr))
}

test_coverage_GRangesList <- function() {
    grl <- make_test_GRangesList()

    target <- RleList(chr1=Rle(1:3, c(4, 1, 5)),
                      chr2=Rle(c(1L, 3L, 5L, 6L, 3L), c(1, 1, 1, 7, 3)),
                      chr3=Rle(0:4, c(6, 1, 1, 1, 1)),
                      chr4=Rle(0:6, c(3, 1, 1, 1, 1, 1, 5)),
                      chr5=Rle(0:4, c(9, 1, 1, 1, 1)),
                      compress=FALSE)
    checkIdentical(target, coverage(grl))

    width <- c(chr1=10, chr2=20, chr3=30, chr4=40, chr5=50)
    target <- RleList(chr1=Rle(1:3, c(4, 1, 5)),
                      chr2=Rle(c(1L, 3L, 5L, 6L, 3L, 0L), c(1, 1, 1, 7, 3, 7)),
                      chr3=Rle(c(0:4, 0L), c(6, 1, 1, 1, 1, 20)),
                      chr4=Rle(c(0:6, 0L), c(3, 1, 1, 1, 1, 1, 5, 27)),
                      chr5=Rle(c(0:4, 0L), c(9, 1, 1, 1, 1, 37)),
                      compress=FALSE)
    checkIdentical(target, coverage(grl, width=width))

    weight <- list(chr1=1L, chr2=10L, chr3=100L, chr4=1000L, chr5=10000L)
    target <- RleList(chr1=Rle(1:3, c(4, 1, 5)),
                      chr2=Rle(10L * c(1L, 3L, 5L, 6L, 3L), c(1, 1, 1, 7, 3)),
                      chr3=Rle(100L * 0:4, c(6, 1, 1, 1, 1)),
                      chr4=Rle(1000L * 0:6, c(3, 1, 1, 1, 1, 1, 5)),
                      chr5=Rle(10000L * 0:4, c(9, 1, 1, 1, 1)),
                      compress=FALSE)
    checkIdentical(target, coverage(grl, weight=weight))

    shift <- list(chr1=0, chr2=1, chr3=2, chr4=3, chr5=4)
    target <- RleList(chr1=Rle(1:3, c(4, 1, 5)),
                      chr2=Rle(c(0L, 1L, 3L, 5L, 6L, 3L), c(1, 1, 1, 1, 7, 3)),
                      chr3=Rle(0:4, c(8, 1, 1, 1, 1)),
                      chr4=Rle(0:6, c(6, 1, 1, 1, 1, 1, 5)),
                      chr5=Rle(0:4, c(13, 1, 1, 1, 1)),
                      compress=FALSE)
    checkIdentical(target, coverage(grl, shift=shift))
}

