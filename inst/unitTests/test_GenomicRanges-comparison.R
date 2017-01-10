###

make_test_GRanges <- function()
    GRanges("chr1",
            IRanges(c(11:13, 13:10, 11:12), width=5, names=LETTERS[1:9]),
            Rle(strand(c("+", "-", "+", "-")), c(4, 3, 1, 1)))

test_is.unsorted_GenomicRanges <- function()
{
    gr <- make_test_GRanges()

    ## Warning expected
    checkException(S4Vectors:::errorIfWarning(is.unsorted(gr, na.rm=TRUE)),
                   silent=TRUE)

    checkTrue(is.unsorted(gr))
    checkTrue(is.unsorted(gr, strictly=TRUE))

    sorted_gr <- sort(gr)
    checkTrue(!is.unsorted(sorted_gr))
    checkTrue(is.unsorted(sorted_gr, strictly=TRUE))
    checkTrue(!is.unsorted(unique(sorted_gr), strictly=TRUE))

    ## Ignore the strand
    checkTrue(is.unsorted(gr, ignore.strand=TRUE))
    checkTrue(is.unsorted(gr, strictly=TRUE, ignore.strand=TRUE))

    sorted_gr <- sort(gr, ignore.strand=TRUE)
    checkTrue(is.unsorted(sorted_gr))
    checkTrue(!is.unsorted(sorted_gr, ignore.strand=TRUE))
    checkTrue(is.unsorted(sorted_gr, strictly=TRUE, ignore.strand=TRUE))

    gr2 <- sorted_gr[c(1:2, 7:8)]
    checkTrue(!is.unsorted(gr2, strictly=TRUE, ignore.strand=TRUE))
    checkTrue(is.unsorted(gr2, strictly=TRUE))
}

test_order_GenomicRanges <- function()
{
    gr <- make_test_GRanges()

    target <- c(1L, 8L, 2L, 3L, 4L, 7L, 6L, 5L, 9L)
    checkTrue(!is.unsorted(gr[target]))
    checkIdentical(target, order(gr))

    target <- c(5L, 9L, 6L, 7L, 3L, 4L, 2L, 1L, 8L)
    checkTrue(!is.unsorted(gr[rev(target)]))
    checkIdentical(target, order(gr, decreasing=TRUE))
}

test_sort_GenomicRanges <- function()
{
    gr <- make_test_GRanges()

    sorted_names <- c("A", "H", "B", "C", "D", "G", "F", "E", "I")
    checkIdentical(gr[sorted_names], sort(gr))

    sorted_names <- c("E", "I", "F", "G", "C", "D", "B", "A", "H")
    checkIdentical(gr[sorted_names], sort(gr, decreasing=TRUE))

    ## Ignore the strand
    sorted_names <- names(sort(unstrand(gr)))
    checkIdentical(gr[sorted_names], sort(gr, ignore.strand=TRUE))

    sorted_names <- names(sort(unstrand(gr), decreasing=TRUE))
    checkIdentical(gr[sorted_names], sort(gr, decreasing=TRUE,
                                              ignore.strand=TRUE))
}

test_rank_GenomicRanges <- function()
{
    gr <- make_test_GRanges()

    target <- c(1L, 3L, 4L, 4L, 8L, 7L, 6L, 1L, 8L)
    checkIdentical(target, rank(gr, ties.method="min"))

    checkIdentical(rank(target), rank(gr))
    checkIdentical(rank(target, ties.method="average"),
                   rank(gr, ties.method="average"))
    checkIdentical(rank(target, ties.method="first"),
                   rank(gr, ties.method="first"))
    checkIdentical(rank(target, ties.method="last"),
                   rank(gr, ties.method="last"))
    checkIdentical(rank(target, ties.method="max"),
                   rank(gr, ties.method="max"))
    checkIdentical(rank(target, ties.method="min"),
                   rank(gr, ties.method="min"))

    ## Ignore the strand
    checkIdentical(rank(unstrand(gr)),
                   rank(gr, ignore.strand=TRUE))
    checkIdentical(rank(unstrand(gr), ties.method="average"),
                   rank(gr, ties.method="average", ignore.strand=TRUE))
    checkIdentical(rank(unstrand(gr), ties.method="first"),
                   rank(gr, ties.method="first", ignore.strand=TRUE))
    checkIdentical(rank(unstrand(gr), ties.method="last"),
                   rank(gr, ties.method="last", ignore.strand=TRUE))
    checkIdentical(rank(unstrand(gr), ties.method="max"),
                   rank(gr, ties.method="max", ignore.strand=TRUE))
    checkIdentical(rank(unstrand(gr), ties.method="min"),
                   rank(gr, ties.method="min", ignore.strand=TRUE))
}

