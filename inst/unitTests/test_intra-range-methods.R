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

test_shift_GenomicRanges <- function()
{
    ## empty, reversibility, recycling 'x'
    gr <- make_test_GRanges()
    checkIdentical(shift(GRanges(), 10), GRanges())
    checkIdentical(gr, shift(shift(gr, 10), -10))
    x <- 1:2
    checkIdentical(start(shift(gr[1:4], x)), start(gr[1:4]) + x)

    ## no seqlength or circularity
    checkIdentical(start(gr) + 10L, start(shift(gr, 10)))
    checkIdentical(width(gr), width(shift(gr, 10)))
    gr <- GRanges("chrA", IRanges(20, 30))
    checkIdentical(IRanges(8, 18), ranges(shift(gr, -12)))
    checkIdentical(IRanges(98, 108), ranges(shift(gr, 78)))

    ## seqlength and circularity combos
    gr <- GRanges("chr1", IRanges(5, width=6))
    isCircular(gr) <- TRUE
    checkIdentical(start(shift(gr, -10)), -5L)

    seqlengths(gr) <- 20
    isCircular(gr) <- NA
    warn <- FALSE
    res <- withCallingHandlers({
        shift(gr, -10)
    }, warning=function(w) {
        warn <<- TRUE
        invokeRestart("muffleWarning")
    })
    checkTrue(warn == TRUE)
    checkIdentical(start(res), -5L)

    isCircular(gr) <- FALSE
    warn <- FALSE
    res <- withCallingHandlers({
        shift(gr, -10)
    }, warning=function(w) {
        warn <<- TRUE
        invokeRestart("muffleWarning")
    })
    checkTrue(warn == TRUE)
    checkIdentical(start(res), -5L)
}

test_shift_GRangesList <- function()
{
    grl <- make_test_GRangesList()
    shifted <- shift(grl, 10)
    checkIdentical(start(grl) + 10L, start(shifted))
}

test_resize_GenomicRanges <- function()
{
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
    ## No-ops.
    checkIdentical(gr, resize(gr, width=width(gr)))
    checkIdentical(gr, resize(gr, width=width(gr), fix="end"))
    checkIdentical(gr, resize(gr, width=width(gr), fix="center"))
}

test_resize_GRangesList <- function()
{
    grl <- make_test_GRangesList()
    target <- endoapply(grl, resize, width=5)
    current <- resize(grl, width=5)
    checkIdentical(target, current)

    ## No-ops.
    checkIdentical(grl, resize(grl, width=width(grl)))
    checkIdentical(grl, resize(grl, width=width(grl), fix="end"))
    checkIdentical(grl, resize(grl, width=width(grl), fix="center"))
}

test_flank_GenomicRanges <- function()
{
    checkIdentical(flank(GRanges(), 10), GRanges())

    gr_seqnames <- c("chr1", "chr2", "chr1", "chrM")
    gr_ranges <- IRanges(21:24, width=10)
    gr_strand <- strand(c("+", "-", "*", "-"))
    gr <- GRanges(gr_seqnames, gr_ranges, gr_strand)

    ## NO warning expected.
    S4Vectors:::errorIfWarning(current <- flank(gr, 10))
    checkTrue(S4Vectors:::errorIfWarning(validObject(current)))
    target_ranges <- IRanges(c(11, 32, 13, 34), width=10)
    target <- GRanges(gr_seqnames, target_ranges, gr_strand)
    checkIdentical(target, current)

    ## NO warning expected.
    S4Vectors:::errorIfWarning(current <- flank(gr, 10, start=FALSE))
    checkTrue(S4Vectors:::errorIfWarning(validObject(current)))
    target_ranges <- IRanges(c(31, 12, 33, 14), width=10)
    target <- GRanges(gr_seqnames, target_ranges, gr_strand)
    checkIdentical(target, current)

    ## NO warning expected.
    S4Vectors:::errorIfWarning(current <- flank(gr, 30))
    checkTrue(S4Vectors:::errorIfWarning(validObject(current)))
    target_ranges <- IRanges(c(-9, 32, -7, 34), width=30)
    target <- GRanges(gr_seqnames, target_ranges, gr_strand)
    checkIdentical(target, current)

    ## NO warning expected.
    S4Vectors:::errorIfWarning(current <- flank(gr, 30, start=FALSE))
    checkTrue(S4Vectors:::errorIfWarning(validObject(current)))
    target_ranges <- IRanges(c(31, -8, 33, -6), width=30)
    target <- GRanges(gr_seqnames, target_ranges, gr_strand)
    checkIdentical(target, current)

    seqlengths(gr) <- c(chr1=60, chr2=50, chrM=35)

    ## Warning expected.
    checkException(S4Vectors:::errorIfWarning(
                       current <- flank(gr, 10)
                   ), silent=TRUE)
    suppressWarnings(current <- flank(gr, 10))

    checkException(S4Vectors:::errorIfWarning(
                       validObject(current)
                   ), silent=TRUE)
    checkTrue(suppressWarnings(validObject(current)))

    target_ranges <- IRanges(c(11, 32, 13, 34), width=10)
    checkIdentical(target_ranges, ranges(current))

    isCircular(gr) <- c(chr1=NA, chr2=FALSE, chrM=TRUE)

    ## NO warning expected.
    S4Vectors:::errorIfWarning(current <- flank(gr, 10))
    checkTrue(S4Vectors:::errorIfWarning(validObject(current)))
    target_ranges <- IRanges(c(11, 32, 13, 34), width=10)
    checkIdentical(target_ranges, ranges(current))

    ## Warning expected.
    checkException(S4Vectors:::errorIfWarning(
                       current <- flank(gr, 20)
                   ), silent=TRUE)
    suppressWarnings(current <- flank(gr, 20))

    checkException(S4Vectors:::errorIfWarning(
                       validObject(current)
                   ), silent=TRUE)
    checkTrue(suppressWarnings(validObject(current)))

    target_ranges <- IRanges(c(1, 32, 3, 34), width=20)
    checkIdentical(target_ranges, ranges(current))
}

test_promoters_GenomicRanges <- function()
{
    checkTrue(length(promoters(GRanges())) == 0)

    ## upstream / downstream
    gr <- GRanges("chr1", IRanges(c(5, 10), width=1), "+")
    target <- GRanges("chr1", IRanges(c(5, 10), width=0), "+")
    current <- promoters(gr, 0, 0)
    checkIdentical(target, current)
    strand(gr) <- c("+", "-")
    target <- IRanges(c(3, 11), width=2)
    current <- ranges(promoters(gr, 2, 0))
    checkIdentical(target, current)
    target <- IRanges(c(5, 9), width=2)
    current <- ranges(promoters(gr, 0, 2))
    checkIdentical(target, current)

    gr <- GRanges("chr1", IRanges(0, width=6), "+")
    target <- GRanges("chr1", IRanges(-3, 2), "+")
    current <- promoters(gr, 3, 3)
    checkIdentical(target, current)
    checkTrue(validObject(current) == TRUE)
    gr <- GRanges("chr1", IRanges(rep(10, 3), width=6), c("+", "-", "*"))
    target <- GRanges("chr1", IRanges(c(7, 13, 7), c(12, 18, 12)),
        c("+", "-", "*"))
    current <- suppressWarnings(promoters(gr, 3, 3))
    checkIdentical(target, current)

    ## treat "*" as "+"
    gr <- GRanges("chr1", IRanges(5, width=6), "+")
    target <- GRanges("chr1", IRanges(2, 7), "+")
    current <- promoters(gr, 3, 3)
    checkIdentical(target, current)
    strand(gr) <- "*"
    strand(target) <- "*"
    current <- suppressWarnings(promoters(gr, 3, 3))
    checkIdentical(target, current)

    ## metadata
    gr <- GRanges("chr1", IRanges(0, width=6), names="A", strand="+", score=99)
    current <- promoters(gr, 3, 3)
    checkIdentical(mcols(gr), mcols(current))
    checkIdentical(names(gr), names(current))
    checkIdentical(seqinfo(gr), seqinfo(current))
}

test_restrict_GenomicRanges <- function()
{
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

test_trim_GenomicRanges <- function()
{
    checkIdentical(trim(GRanges()), GRanges())

    gr_seqnames <- c("chr1", "chr2", "chr1", "chrM")
    gr_ranges <- IRanges(0:3, width=30)

    ## NO warning expected.
    S4Vectors:::errorIfWarning(gr <- GRanges(gr_seqnames, gr_ranges))
    checkTrue(S4Vectors:::errorIfWarning(validObject(gr)))
    checkIdentical(trim(gr), gr)

    gr_seqlengths <- c(chr1=50, chr2=NA, chrM=NA)

    ## Warning expected.
    checkException(S4Vectors:::errorIfWarning(
                       seqlengths(gr) <- gr_seqlengths
                   ), silent=TRUE)
    suppressWarnings(seqlengths(gr) <- gr_seqlengths)

    checkException(S4Vectors:::errorIfWarning(
                       validObject(gr)
                   ), silent=TRUE)
    checkTrue(suppressWarnings(validObject(gr)))

    gr <- trim(gr)
    checkTrue(S4Vectors:::errorIfWarning(validObject(gr)))
    target_ranges <- IRanges(c(1, 1, 2, 3), width=c(29, 30, 30, 30))
    checkIdentical(target_ranges, ranges(gr))

    isCircular(gr) <- c(chr1=FALSE, chr2=FALSE, chrM=TRUE)

    ## NO warning expected.
    gr_seqlengths <- c(chr1=50, chr2=NA, chrM=15)
    S4Vectors:::errorIfWarning(seqlengths(gr) <- gr_seqlengths)
    checkTrue(S4Vectors:::errorIfWarning(validObject(gr)))
    checkIdentical(trim(gr), gr)
}

