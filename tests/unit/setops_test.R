make_test_GRanges <- function() {
    new("GRanges",
        seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        values = DataFrame(score = 1:10, GC = seq(1, 0, length=10)))
}

make_test_GRangesList <- function() {
    GRangesList(
        a =
        new("GRanges",
            seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
            ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
            strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
            values = DataFrame(score = 1:10, GC = seq(1, 0, length=10))),
        b =
        new("GRanges",
            seqnames = Rle(c("chr2", "chr4", "chr5"), c(3, 6, 4)),
            ranges = IRanges(1:13, width = 13:1, names = tail(letters, 13)),
            strand = Rle(strand(c("-", "+", "-")), c(4, 5, 4)),
            values = DataFrame(score = 1:13, GC = seq(0, 1, length=13))))
}

test_pintersect <- function()
{
    gr <- make_test_GRanges()
    grlist <- make_test_GRangesList()

    checkException(pintersect(grlist, grlist), silent = TRUE)

    expect <- gr
    values(expect) <- NULL
    checkIdentical(pintersect(gr, gr), expect)
    checkIdentical(pintersect(gr[1:2], grlist), pintersect(grlist, gr[1:2]))

    expect <-
      GRangesList(a =
                  GRanges(c("chr2", "chr2", "chr2"),
                          IRanges(c(4,4,4), c(10,10,10), names = letters[2:4]),
                          c("+", "+", "*")),
                  b = GRanges())
    checkIdentical(pintersect(gr[4:5], grlist), expect)
}
