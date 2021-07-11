make_test_GRanges <- function()
    GRanges(Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
            IRanges(1:10, width=10:1, names=head(letters, 10)),
            Rle(c("-", "+", "*", "+", "-"), c(1, 2, 2, 3, 2)),
            score=1:10, GC=seq(1, 0, length=10),
            seqinfo=Seqinfo(paste("chr", 1:3, sep="")))

test_range_GenomicRanges <- function()
{
    gr <- make_test_GRanges()
    current <- range(gr)
    target <- GRanges(Rle(c("chr1", "chr2", "chr3"), c(3, 2, 2)),
                      IRanges(start=c(6, 1, 5, 2, 4, 7, 9), end=10),
                      c("+", "-", "*", "+", "*", "+", "-"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, current)

    current <- range(gr, ignore.strand=TRUE)
    target <- GRanges(c("chr1", "chr2", "chr3"),
                      IRanges(start=c(1, 2, 7), end=10),
                      c("*", "*", "*"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, current)

    # test with.revmap
    current <- range(gr, with.revmap=TRUE, ignore.strand=TRUE)
    mcols(target)$revmap <- IntegerList(c(1,5,6),  c(2,3,4),  c(7:10))
    checkIdentical(target, current)
}

test_range_GRangesList <- function()
{
    gr <- make_test_GRanges()
    grl <- GRangesList(gr, shift(rev(gr), 5 * seq_along(gr)))
    for (ignore.strand in c(FALSE, TRUE)) {
        current <- range(grl, ignore.strand=TRUE)
        target <- endoapply(grl, range, ignore.strand=TRUE)
        checkTrue(validObject(current, complete=TRUE))
        checkIdentical(target, current)
    }

    # test with.revmap
    obj <- range(grl, with.revmap=TRUE, ignore.strand=TRUE)
    revmap1 <- mcols(obj[[1]])$revmap
    revmap2 <- mcols(obj[[2]])$revmap
    ans1 <- IntegerList(c(1,5,6),  c(2,3,4),  c(7:10))
    ans2 <- IntegerList(c(5,6,10), c(7:9), c(1:4))
    checkIdentical(revmap1, ans1)
    checkIdentical(revmap2, ans2)
}

test_reduce_GenomicRanges <- function()
{
    gr <- make_test_GRanges()
    current <- reduce(gr)
    target <- GRanges(Rle(c("chr1", "chr2", "chr3"), c(3, 2, 2)),
                      IRanges(start=c(6, 1, 5, 2, 4, 7, 9), end=10),
                      c("+", "-", "*", "+", "*", "+", "-"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, current)

    current <- reduce(gr, with.revmap=TRUE)
    mcols(target)$revmap <- IntegerList(6, 1, 5, 2:3, 4, 7:8, 9:10)
    checkIdentical(target, current)
}

test_reduce_GRangesList <- function()
{
    gr <- make_test_GRanges()
    grl <- GRangesList(gr, shift(rev(gr), 5 * seq_along(gr)))
    for (with.revmap in c(FALSE, TRUE)) {
        for (ignore.strand in c(FALSE, TRUE)) {
            current <- reduce(grl, with.revmap=with.revmap,
                                   ignore.strand=ignore.strand)
            target <- endoapply(grl, reduce, with.revmap=with.revmap,
                                             ignore.strand=ignore.strand)
            checkTrue(validObject(current, complete=TRUE))
            checkIdentical(target, current)
        }
    }
}

test_gaps_GenomicRanges <- function()
{
    gr <- make_test_GRanges()
    current <- gaps(gr, start=1, end=10)
    target <- GRanges(Rle(c("chr1", "chr2", "chr3"), c(2, 3, 3)),
                      IRanges(start=1, end=c(5, 4, 1, 10, 3, 6, 8, 10)),
                      c("+", "*", "+", "-", "*", "+", "-", "*"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, current)
}

test_disjoin_GenomicRanges <- function()
{
    gr <- make_test_GRanges()
    current <- disjoin(gr)
    target <- GRanges(Rle(c("chr1", "chr2", "chr3"), c(3, 3, 4)),
                      IRanges(start=c(6, 1, 5, 2, 3, 4, 7, 8, 9, 10),
                              end=c(10, 10, 10, 2, 10, 10, 7, 10, 9, 10)),
                      c("+", "-", "*", "+", "+", "*", "+", "+", "-", "-"))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, current)

    gr <- GRanges(Rle(c("chr1", "chr3"), c(2, 2)),
                  IRanges(c(8, 6, 8, 6), c(11, 15, 11, 15),
                          names=c("k", "l", "m", "n")),
                  c("-", "-", "+", "*"),
                  score=11:14, GC=c(.2, .3, .3, .1))
    current <- disjoin(gr)
    target <- GRanges(Rle(c("chr1", "chr3"), c(3, 2)),
                      IRanges(c(6, 8, 12, 8, 6), c(7, 11, 15, 11, 15)),
                      Rle(c("-", "+", "*"), c(3, 1, 1)))
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical(target, current)

    current <- disjoin(gr, with.revmap=TRUE)
    mcols(target)$revmap <- IntegerList(2, 1:2, 2, 3, 4)
    checkIdentical(target, current)
}

test_disjoin_GRangesList <- function()
{
    grl <- GRangesList(make_test_GRanges(),
                       GRanges("1", IRanges(1, 10), score=21, GC=.21),
                       GRanges(),
                       GRanges(Rle(c("chr1", "chr3"), c(2, 2)),
                               IRanges(c(8, 6, 8, 6), c(11, 15, 11, 15),
                                       names=c("k", "l", "m", "n")),
                               strand(c("-", "-","+","*")),
                               score=41:44, GC=c(.41, .42, .43, .44)))
    for (with.revmap in c(FALSE, TRUE)) {
        for (ignore.strand in c(FALSE, TRUE)) {
            current <- disjoin(grl, with.revmap=with.revmap,
                                    ignore.strand=ignore.strand)
            target <- endoapply(grl, disjoin, with.revmap=with.revmap,
                                              ignore.strand=ignore.strand)
            checkTrue(validObject(current, complete=TRUE))
            checkIdentical(target, current)
        }
    }
}

