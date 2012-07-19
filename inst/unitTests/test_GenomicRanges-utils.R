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
    checkIdentical(elementMetadata(gr), elementMetadata(res))
    checkIdentical(seqnames(gr), seqnames(res))
    checkIdentical(seqinfo(gr), seqinfo(res))
    target <- IRanges(start=c(4, 5, 5, 5, 5, 6, 7, 8, 9, 10),
                      end = c(10, 8, 8, 8, 10, 10, 9, 9, 9, 9),
                      names=letters[1:10])
    checkIdentical(ranges(res), target)
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

test_GenomicRanges_findNearest0 <- function()
{
    .findNearest <- GenomicRanges:::.GenomicRanges_findNearest0
    sentinel <- c(0, 20) 
    subject <- c(5, 15) 

    hits <- .findNearest(10, subject, sentinel, TRUE)
    checkIdentical(2L, subjectHits(hits))
    hits <- .findNearest(10, subject, sentinel, FALSE)
    checkIdentical(1L, subjectHits(hits))

    hits <- .findNearest(5, subject, sentinel, TRUE)
    checkIdentical(2L, subjectHits(hits))
    hits <- .findNearest(15, subject, sentinel, FALSE)
    checkIdentical(1L, subjectHits(hits))
 
    hits <- .findNearest(5, subject, sentinel, FALSE)
    checkIdentical(integer(), subjectHits(hits))
    hits <- .findNearest(15, subject, sentinel, TRUE)
    checkIdentical(integer(), subjectHits(hits))

    subject <- c(15, 5)
    hits <- .findNearest(10, subject, sentinel, TRUE)
    checkIdentical(1L, subjectHits(hits))
    hits <- .findNearest(10, subject, sentinel, FALSE)
    checkIdentical(2L, subjectHits(hits))
}

test_GenomicRanges_precede_follow <- function()
{
    ## query on "+"
    query <- GRanges("A", IRanges(c(1, 5, 10, 15, 20), width=1), "+")
    subject <- GRanges("A", IRanges(c(5, 15), width=1), "+")
    hits <- precede(query, subject)
    checkIdentical(c(1L, 2L, 2L, NA_integer_, NA_integer_), hits)
    hits <- follow(query, subject)
    checkIdentical(c(NA_integer_, NA_integer_, 1L, 1L, 2L), hits)

    subject <- GRanges("A", IRanges(c(5, 15), width=1), "-")
    hits <- precede(query, subject)
    checkIdentical(rep(NA_integer_, length(query)), hits)
    hits <- follow(query, subject)
    checkIdentical(rep(NA_integer_, length(query)), hits)

    subject <- GRanges("A", IRanges(c(5, 15), width=1), "*")
    hits <- precede(query, subject)
    checkIdentical(c(1L, 2L, 2L, NA_integer_, NA_integer_), hits)
    hits <- follow(query, subject)
    checkIdentical(c(NA_integer_, NA_integer_, 1L, 1L, 2L), hits)

    ## query on "-"
    query <- GRanges("A", IRanges(c(1, 5, 10, 15, 20), width=1), "-")
    subject <- GRanges("A", IRanges(c(5, 15), width=1), "-")
    hits <- precede(query, subject)
    checkIdentical(c(NA_integer_, NA_integer_, 1L, 1L, 2L), hits)
    hits <- follow(query, subject)
    checkIdentical(c(1L, 2L, 2L, NA_integer_, NA_integer_), hits)

    subject <- GRanges("A", IRanges(c(5, 15), width=1), "+")
    hits <- precede(query, subject)
    checkIdentical(rep(NA_integer_, length(query)), hits)
    hits <- follow(query, subject)
    checkIdentical(rep(NA_integer_, length(query)), hits)

    subject <- GRanges("A", IRanges(c(5, 15), width=1), "*")
    hits <- precede(query, subject)
    checkIdentical(c(NA_integer_, NA_integer_, 1L, 1L, 2L), hits)
    hits <- follow(query, subject)
    checkIdentical(c(1L, 2L, 2L, NA_integer_, NA_integer_), hits)

    ## query on "*"
    query <- GRanges("A", IRanges(c(1, 5, 10, 15, 20), width=1), "*")
    subject <- GRanges("A", IRanges(c(5, 15), width=1), "+")
    hits <- precede(query, subject)
    checkIdentical(c(1L, 2L, 2L, NA_integer_, NA_integer_), hits)
    hits <- follow(query, subject)
    checkIdentical(c(NA_integer_, NA_integer_, 1L, 1L, 2L), hits)

    subject <- GRanges("A", IRanges(c(5, 15), width=1), "-")
    hits <- precede(query, subject)
    checkIdentical(c(NA_integer_, NA_integer_, 1L, 1L, 2L), hits)
    hits <- follow(query, subject)
    checkIdentical(c(1L, 2L, 2L, NA_integer_, NA_integer_), hits)

    subject <- GRanges("A", IRanges(c(5, 15), width=1), "*")
    hits <- precede(query, subject)
    checkIdentical(c(1L, 2L, 1L, 1L, 2L), hits)
    hits <- follow(query, subject)
    checkIdentical(c(1L, 2L, 1L, 1L, 2L), hits)
}

test_GenomicRanges_precede_follow_ties <- function()
{
    query <- GRanges("A", IRanges(10, width=1), c("+", "-", "*"))
    subject <- GRanges("A", IRanges(c(5, 5, 5, 15, 15, 15), width=1),
                       rep(c("+", "-", "*"), 2))
    checkIdentical(c(4L, 2L, 2L), precede(query, subject))
    checkIdentical(c(1L, 4L, 1L), precede(query, rev(subject)))
}

test_GenomicRanges_ignore_strand <- function()
{
    query <- GRanges("A", IRanges(10, width=1), c("+", "-", "*"))
    subject <- GRanges("A", IRanges(c(5, 5, 5, 15, 15, 15), width=1),
                       rep(c("+", "-", "*"), 2))
    checkIdentical(c(4L, 4L, 4L),
                   precede(query, subject, ignore.strand=TRUE))
    checkIdentical(c(1L, 1L, 1L), 
                   precede(query, rev(subject),ignore.strand=TRUE))
}

test_GenomicRanges_nearest <- function()
{
    r <- IRanges(c(1,5,10), c(2,7,12))
    g <- GRanges("chr1", r, "+")
    checkEquals(precede(r), precede(g))
    checkEquals(follow(r), follow(g))
    checkEquals(nearest(r), nearest(g))

    g <- GRanges("chr1", r, "-")
    checkEquals(follow(r), precede(g))
    checkEquals(precede(r), follow(g))
    checkEquals(nearest(r), nearest(g))

    g <- GRanges("chr1", r, "*")
    checkEquals(follow(g), precede(g))
    checkEquals(nearest(r), follow(g))
    checkEquals(follow(g), nearest(g))
}

test_GenomicRanges_distance <- function()
{
    g1 <- GRanges(seqnames = c(rep("chr1", 3), rep("chr2", 2)),
        ranges = IRanges(rep(1, 5),  width=3),
        strand = c("+", "-", "*", "*", "*"))

    g2 <- GRanges(seqnames = c(rep("chr1", 3), rep("chr2", 2)),
        ranges = IRanges(rep(5, 5),  width=3),
        strand = c("+", "-", "*", "-", "+"))

    current <- distance(g1, g2)
    target <- c(2L, 2L, 2L, 2L, 2L)
    checkIdentical(current, target)

    strand(g2[1]) <- "-"
    current <- distance(g1[1], g2[1])
    checkTrue(is.na(current))

    seqnames(g2[4]) <- factor("chr1", levels=c("chr1", "chr2")) 
    current <- distance(g1[4], g2[4])
    checkTrue(is.na(current))
}

