###

make_test_GRanges <- function() {
    ## TODO: Make a GRanges object that is more challenging to check whether
    ## it is sorted.
    GRanges(
        Rle(factor(c("chr1", "chr2", "chr1", "chr3"),
                   levels=paste0("chr", 1:5)),
            c(1, 3, 2, 4)),
        IRanges(1:10, width=10:1, names=head(letters, 10)),
        strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        score=1:10,
        GC=seq(1, 0, length=10)
   )
}

test_is.unsorted <- function()
{
    gr <- make_test_GRanges()

    ## Warning expected
    checkException(S4Vectors:::errorIfWarning(is.unsorted(gr, na.rm=TRUE)),
                   silent=TRUE)

    checkTrue(is.unsorted(gr))
    checkTrue(!is.unsorted(sort(gr)))

    ## is.unsorted() should return TRUE if 'strictly=TRUE' and 'gr' contains
    ## duplicate ranges, even if these are sorted.
    dup_gr <- sort(c(gr[1], gr[1]))
    checkTrue(!is.unsorted(dup_gr))
    checkTrue(is.unsorted(dup_gr, strictly=TRUE))

    ## Ranges that differ in only seqnames, strand, start, or width
    x <- c(gr[1], gr[1])
    xx <- x
    seqnames(xx) <- factor(c("chr2", "chr1"), levels=seqlevels(xx))
    checkTrue(is.unsorted(xx))
    checkTrue(!is.unsorted(rev(xx)))

    xx <- x
    strand(xx) <- c("-", "+")
    checkTrue(is.unsorted(xx))
    checkTrue(!is.unsorted(rev(xx)))

    xx <- x
    start(xx) <- c(2, 1)
    checkTrue(is.unsorted(xx))
    checkTrue(!is.unsorted(rev(xx)))

    xx <- x
    width(xx) <- c(10, 7)
    checkTrue(is.unsorted(xx))
    checkTrue(!is.unsorted(rev(xx)))

    ## Ignore the strand
    checkTrue(is.unsorted(gr, ignore.strand=TRUE))
    gr2 <- sort(gr, ignore.strand=TRUE)
    checkTrue(is.unsorted(gr2))
    checkTrue(!is.unsorted(gr2, ignore.strand=TRUE))
}

