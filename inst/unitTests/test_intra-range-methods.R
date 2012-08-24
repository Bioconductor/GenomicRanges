make_test_GRanges <- function() {
    new("GRanges",
        seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")),
        elementMetadata = DataFrame(score = 1:10, GC = seq(1, 0, length=10)))
}

test_GenomicRanges_intra_interval_ops <- function()
{
    ## shift
    gr <- make_test_GRanges()
    shifted <- shift(gr, 10)
    checkIdentical(start(gr) + 10L, start(shifted))
    checkIdentical(width(gr), width(shifted))
    gr <- GRanges("chrA", IRanges(20, 30), seqlengths=c(chrA=100))
    checkIdentical(IRanges(8, 18), ranges(shift(gr, -12)))
    shifted <- suppressWarnings(shift(gr, 78))
    checkIdentical(IRanges(98, 100), ranges(shifted))

    ## flank
    gr <- make_test_GRanges()
    flanked <- flank(gr, 10)
    checkIdentical(rep(10L, length(gr)), width(flanked))
    checkIdentical(ifelse(as.vector(strand(gr) != "-"),
                          start(gr) - 10L, end(gr) + 1L), start(flanked))
    flanked <- flank(gr, 10, FALSE)
    checkIdentical(rep(10L, length(gr)), width(flanked))
    checkIdentical(ifelse(as.vector(strand(gr) != "-"),
                          end(gr) + 1L, start(gr) - 10L), start(flanked))

    ## resize
    gr <- make_test_GRanges()
    checkException(resize(gr, 10, fix = "middle"), silent = TRUE)
    checkException(resize(gr, 10, fix = rep("end", 3)), silent = TRUE)
    resized <- resize(gr, 10)
    checkIdentical(rep(10L, length(gr)), width(resized))
    checkIdentical(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 1L), start(resized))
    checkIdentical(ranges(resize(gr, 10, fix = "center")),
                   IRanges(rep(1:5, each=2), width = 10,
                           names = head(letters, 10)))
    checkIdentical(ranges(resize(gr, 10, fix = c("start", "end"))),
                   IRanges(c(1L, 1L, 3L, 1L, 5L, 1L, 7L, 1L, 1L, 10L),
                           width = 10, names = head(letters, 10)))
 
    ## restrict
    gr <-  make_test_GRanges()
    st <- structure(c(4,5), names = c("chr1", "chr2"))
    en <-  structure(c(8,9), names = c("chr2", "chr3"))
    res <- restrict(gr, start = st, end = en)
    checkIdentical(mcols(gr), mcols(res))
    checkIdentical(seqnames(gr), seqnames(res))
    checkIdentical(seqinfo(gr), seqinfo(res))
    target <- IRanges(start=c(4, 5, 5, 5, 5, 6, 7, 8, 9, 10),
                      end = c(10, 8, 8, 8, 10, 10, 9, 9, 9, 9),
                      names=letters[1:10])
    checkIdentical(ranges(res), target)
}

