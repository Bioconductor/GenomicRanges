library(TxDb.Hsapiens.UCSC.hg19.knownGene)

quiet <- suppressWarnings
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

test_GenomicRanges_nearest <- function()
{
    ## adjacent
    r <- IRanges(c(1,6), c(5,10))
    g <- GRanges("chr1", r, "+")
    checkEquals(follow(r), follow(g))
    checkEquals(precede(r), precede(g))
    checkEquals(nearest(r), nearest(g))
    checkEquals(follow(r,r), follow(g,g))
    checkEquals(precede(r,r), precede(g,g))
    checkEquals(nearest(r,r), nearest(g,g))
    g <- GRanges("chr1", r, "-")
    checkEquals(follow(r), precede(g))
    checkEquals(precede(r), follow(g))
    checkEquals(nearest(r), nearest(g))
    checkEquals(follow(r,r), precede(g,g))
    checkEquals(precede(r,r), follow(g,g))
    checkEquals(nearest(r,r), nearest(g,g))
    g <- GRanges("chr1", r, "*")
    checkEquals(nearest(r), nearest(g))
    checkEquals(nearest(r,r), nearest(g,g))

    ## separated by 1
    r <- IRanges(c(1,7), c(5,11))
    g <- GRanges("chr1", r, "+")
    checkEquals(follow(r), follow(g))
    checkEquals(precede(r), precede(g))
    checkEquals(nearest(r), nearest(g))
    checkEquals(nearest(r,r), nearest(g,g))
    g <- GRanges("chr1", r, "-")
    checkEquals(follow(r), precede(g))
    checkEquals(precede(r), follow(g))
    checkEquals(nearest(r), nearest(g))
    checkEquals(nearest(r,r), nearest(g,g))
    g <- GRanges("chr1", r, "*")
    checkEquals(nearest(r), nearest(g))
    checkEquals(nearest(r,r), nearest(g,g))

    ## separated by > 1
    r <- IRanges(c(1,5,10), c(2,7,12))
    g <- GRanges("chr1", r, "+")
    checkEquals(precede(r), precede(g))
    checkEquals(follow(r), follow(g))
    checkEquals(nearest(r), nearest(g))
    g <- GRanges("chr1", r, "-")
    checkEquals(follow(r), precede(g))
    checkEquals(precede(r), follow(g))
    checkEquals(nearest(r), nearest(g))

    ## '*' strand precedes or follows both ranges
    x <- GRanges("chr1", IRanges(1100, width=1), strand='*')
    y <- GRanges("chr1", IRanges(c(1000, 1500), width=1), strand=c("+", "-"))
    checkTrue(nearest(x, y) == 1L)
    y <- GRanges("chr1", IRanges(c(1000, 1500), width=1), strand=c("-", "+"))
    checkTrue(nearest(x, y) == 1L)

    ## overlapping
    r <- IRanges(c(1,4,8), c(6,10,12))
    g <- GRanges("chr1", r, "+")
    checkEquals(nearest(r), nearest(g))
    checkEquals(nearest(r, r), nearest(g, g))
    checkEquals(nearest(r, rev(r)), nearest(g, rev(g)))
    g <- GRanges("chr1", r, "-")
    checkEquals(nearest(r), nearest(g))
    checkEquals(nearest(r, r), nearest(g, g))
checkEquals(nearest(r, rev(r)), nearest(g, rev(g)))
g <- GRanges("chr1", r, "*")
    checkEquals(nearest(r), nearest(g))
    checkEquals(nearest(r, r), nearest(g, g))
    checkEquals(nearest(r, rev(r)), nearest(g, rev(g)))

    q <- GRanges("chr1", IRanges(1, 15), "+")
    s <- GRanges("chr1", IRanges(c(1, 1, 10), c(5, 15, 15)), "+")
    target <- nearest(q, s, select="arbitrary")
    checkEquals(3, target)
    strand(q) <- "-"
    strand(s) <- "-"
    target <- nearest(q, s, select="arbitrary")
    checkEquals(3, target)
    target1 <- nearest(ranges(q), ranges(s), select="all")
    target2 <- nearest(q, s, select="all")
    checkEquals(target1, target2)

    ## select = 'all'
    q <- GRanges("chr1", IRanges(1, 1))
    s <- GRanges("chr1", IRanges(5, 5), strand="-")
    target1 <- nearest(ranges(q), ranges(s), select="all")
    target2 <- nearest(q, s, select="all")
    checkEquals(target1, target2)

    ## ignore.strand
    q <- GRanges("chr1", IRanges(5, width=1), "+")
    s <- GRanges("chr1", IRanges(c(10, 8), width=1), "-")
    res <- nearest(q, s, ignore.strand=FALSE)
    checkEquals(res, NA_integer_)
    res <- nearest(q, s, ignore.strand=TRUE)
    checkEquals(res, 2L)

    q <- GRanges("chr1", IRanges(5, 5))
    s <- GRanges("chr1", IRanges(c(6,7), c(6,7)), c("+", "-"))
    checkEquals(nearest(q, s), 1L)

    q <- GRanges("chr1", IRanges(105, 105), "-")
    s <- GRanges("chr1", IRanges(c(1,120), c(100, 125)), c("-", "-"))
    pos <- nearest(q, s, ignore.strand=TRUE)
    neg <- nearest(q, s, ignore.strand=FALSE)
    checkEquals(pos, neg)
}

test_GenomicRanges_distance <- function()
{
    ## empty
    checkException(quiet(distance(GRanges())), silent=TRUE)
    checkIdentical(quiet(distance(GRanges(), GRanges())), integer())

    g1 <- GRanges(seqnames = c(rep("chr1", 3), rep("chr2", 2)),
        ranges = IRanges(rep(1, 5),  width=3),
        strand = c("+", "-", "*", "*", "*"))

    g2 <- GRanges(seqnames = c(rep("chr1", 3), rep("chr2", 2)),
        ranges = IRanges(rep(5, 5),  width=3),
        strand = c("+", "-", "*", "-", "+"))

    current <- quiet(distance(g1, g2))
    target <- rep(1L, length(current))
    checkIdentical(current, target)

    strand(g2[1]) <- "-"
    current <- quiet(distance(g1[1], g2[1]))
    checkTrue(is.na(current))

    seqnames(g2[4]) <- factor("chr1", levels=c("chr1", "chr2"))
    current <- quiet(distance(g1[4], g2[4]))
    checkTrue(is.na(current))

    ## adjacent, overlap, separated by 1
    query <- GRanges("A", IRanges(c(1, 3, 9), c(2, 7, 10)))
    subject <- GRanges("A", IRanges(c(3, 5, 12), c(3, 6, 12)))
    checkIdentical(quiet(distance(query, subject)), c(0L, 0L, 1L))

    ## recycling
    checkIdentical(quiet(distance(query[1:2], subject)),
                   c(0L, 0L, 9L))

    ## zero-width
    target <- abs(-3:3)
    current <- sapply(-3:3, function(i)
                   quiet(distance(shift(IRanges(4,3), i), IRanges(4,3))))
    checkIdentical(current, target)
    query <- GRanges("A", IRanges(4,3))
    subject <-  GRanges("A", IRanges(3,4))
    current <- quiet(distance(query, subject))
    checkIdentical(current, 0L)
}

test_GenomicRanges_distanceToNearest <- function()
{
    target <- Hits(distance=integer(0), sort.by.query=TRUE)
    current <- quiet(distanceToNearest(GRanges()))
    checkIdentical(current, target)

    x <- GRanges("chr1", IRanges(c(1, 5, 10), width=1))
    subject <- GRanges("chr1", IRanges(c(3, 12), width=1))
    current <- quiet(distanceToNearest(x, subject))
    target <- c(1L, 1L, 2L)
    checkIdentical(target, subjectHits(current))

    ## strand
    strand(x) <- "+"
    strand(subject) <- c("+", "-")
    current <- quiet(distanceToNearest(x, subject))
    target <- c(1L, 1L, 1L)
    checkIdentical(target, subjectHits(current))
    current <- quiet(distanceToNearest(x, subject, ignore.strand=TRUE))
    target <- c(1L, 1L, 2L)
    checkIdentical(target, subjectHits(current))

    ## no self-hits / self-hits
    current <- quiet(distanceToNearest(x))
    target <- c(2L, 1L, 2L)
    checkIdentical(target, subjectHits(current))
    current <- quiet(distanceToNearest(x, x))
    target <- c(1L, 2L, 3L)
    checkIdentical(target, subjectHits(current))

    ## ranges start at 0
    x <- GRanges("chr1", IRanges(0, width=1))
    subject <- GRanges("chr1", IRanges(1, width=1))
    current <- distanceToNearest(x, subject)
    iranges <- distanceToNearest(ranges(x), ranges(subject))
    checkIdentical(subjectHits(iranges), subjectHits(current))
    checkIdentical(mcols(iranges)$distance, mcols(current)$distance)
    current <- distanceToNearest(x, x)
    checkIdentical(1L, subjectHits(current))
    checkIdentical(0L, mcols(current)$distance)
    current <- distanceToNearest(subject, x)
    iranges <- distanceToNearest(ranges(subject), ranges(x))
    checkIdentical(subjectHits(iranges), subjectHits(current))
    checkIdentical(mcols(iranges)$distance, mcols(current)$distance)
    current <- distanceToNearest(x, subject)
    iranges <- distanceToNearest(ranges(x), ranges(subject))
    checkIdentical(subjectHits(iranges), subjectHits(current))
    checkIdentical(mcols(iranges)$distance, mcols(current)$distance)

    ## NA handling
    x <- GRanges(c("chr1", "chr2"), IRanges(c(1, 5), width=1))
    current <- quiet(distanceToNearest(x, x[1]))
    checkIdentical(1L, queryHits(current))
    checkIdentical(1L, subjectHits(current))
    checkIdentical(0L, mcols(current)$distance)

    ## select = 'all' vs 'arbitrary'
    x <- GRanges("chr1", IRanges(2, width=1))
    subject <- GRanges("chr1", IRanges(c(8, 10, 8), width=1))
    current <- distanceToNearest(x, subject, select="all")
    checkTrue(is(current, "Hits"))
    checkIdentical(length(current), 2L)
    current <- distanceToNearest(x, GRanges(), select="all")
    checkIdentical(length(current), 0L)
    checkIdentical(queryLength(current), 1L)
}

test_GenomicRanges_findKNN <- function()
{
    ## single argument
    g1 <- GRanges()
    checkIdentical(IntegerList(), findKNN(g1))

    ## empty ranges
    g1 <- GRanges()
    g2 <- GRanges()
    checkIdentical(IntegerList(), findKNN(g1, g2))

    ## no nearest neighbors
    seqinfo <- Seqinfo(paste0("chr", 1:2))
    g1 <- GRanges(c("chr1:1-5:+", "chr1:6-10:+"), seqinfo = seqinfo)
    g2 <- GRanges(c("chr2:1-5:+", "chr2:6-10:+"), seqinfo = seqinfo)
    target <- nearest(g1, g2)
    current <- findKNN(g1, g2)
    checkIdentical(length(target), length(current))
    checkIdentical(target, unlist(current))
    ## NOTE: absolete with now that integer(0) is not NA
    #checkIdentical(lengths(current), c(0L, 0L))

    ## Same strand non-overlapping, x only
    g <- GRanges(c("chr1:11-15:+", "chr1:21-25:+"))
    checkIdentical(nearest(g), unlist(findKNN(g)))

    ## Same strand non-overlapping, x and subject
    g1 <- GRanges(c("chr1:1-5:+", "chr1:7-11:+"))
    g2 <- GRanges(c("chr1:2-3:+", "chr1:10-12:+"))
    checkIdentical(nearest(g1, g2), unlist(findKNN(g1, g2)))

    ## Same strand non-overlapping, some have no neighbors
    g1 <- GRanges(c("chr1:1-5:+", "chr1:11-15:+", "chr2:11-15:+"))
    g2 <- GRanges("chr1:8-9:+")
    checkIdentical(IntegerList(1, 1, NA), findKNN(g1, g2))

    ## Same strand, self, overlapping
    r <- IRanges(c(1,4,8), c(6,10,12))
    g <- GRanges("chr1", r, "+")
    target <- nearest(g)
    near <- findKNN(g)
    ## NOTE: Incorrect due to arbitrarly choosing an overlap
    #checkIdentical(nearest(g), unlist(findKNN(g)))
    #checkIdentical(length(near), 2L)
    #checkIdentical(length(near[[1]]), 1L)
    #checkIdentical(length(near[[2]]), 1L)
    #checkIdentical(near[[1]], 1L)
    #checkIdentical(near[[2]], 2L)
 
    ## Same strand, with overlaps in ranges
    g1 <- GRanges(c("chr1:1-5:+", "chr1:7-11:+"))
    g2 <- GRanges(c("chr1:1-2:+", "chr1:5-7:+", "chr1:10-12:+"))
    target <- nearest(g1, g2)
    near <- findKNN(g1, g2)
    ## NOTE: Incorrect due to arbitrarly choosing an overlap
    #checkIdentical(nearest(g), unlist(findKNN(g)))
    checkIdentical(length(near), 2L)
    checkIdentical(length(near[[1]]), 1L)
    checkIdentical(length(near[[2]]), 1L)
    checkIdentical(near[[1]], 1L)
    checkIdentical(near[[2]], 2L)

    ## Same strand, k > 1
    near <- findKNN(g1, g2, k=3)
    checkIdentical(length(near), 2L)
    checkIdentical(length(near[[1]]), 3L)
    checkIdentical(length(near[[2]]), 3L)
    checkIdentical(near[[1]], c(1L, 2L, 3L))
    checkIdentical(near[[2]], c(2L, 3L, 1L))

    near_all <- findKNN(g1, g2, select="all")
    checkIdentical(length(near_all), 2L)
    checkIdentical(length(near_all[[1]]), 1L)
    checkIdentical(length(near_all[[2]]), 1L)
    checkIdentical(near_all[[1]], c(1L, 2L))
    checkIdentical(near_all[[2]], c(2L, 3L))
    near <- findKNN(g1, g2, k=2)
    checkIdentical(near, near_all)

    ## alternate strands, with overlaps in ranges
    g3 <- GRanges(c("chr1:1-5:-", "chr1:7-11:+"))
    #g4 <- GRanges(c("chr1:1-2:+", "chr1:5-7:-", "chr1:10-12:+"))
    g4 <- GRanges(c("chr1:1-2:-", "chr1:5-7:+", "chr1:10-12:+"))
    near <- findKNN(g3, g4)
    checkIdentical(length(near), 2L)
    checkIdentical(length(near[[1]]), 1L)
    checkIdentical(length(near[[2]]), 1L)
    checkIdentical(near[[1]], 1L)
    checkIdentical(near[[2]], 2L)

    near <- findKNN(g3, g4, k=2)
    checkIdentical(length(near), 2L)
    checkIdentical(length(near[[1]]), 1L)
    checkIdentical(length(near[[2]]), 2L)
    checkIdentical(near[[1]], 1L)
    checkIdentical(near[[2]], c(2L, 3L))

    nearest_all <- findKNN(g3, g4, select="all")
    checkIdentical(length(near_all), 2L)
    checkIdentical(length(near_all[[1]]), 1L)
    checkIdentical(length(near_all[[2]]), 2L)
    checkIdentical(near_all[[1]], 1L)
    checkIdentical(near_all[[2]], c(2L, 3L))
    checkIdentical(near, near_all)

    ## FIXME: Neither nearest() or findKNN() have ignore overlap
    ## ignore.overlap, with overlaps in ranges
#    near <- findKNN(g1, g2, ignore.overlap = TRUE)
#    checkIdentical(length(near), 2L)
#    checkIdentical(length(near[[1]]), 1L)
#    checkIdentical(length(near[[2]]), 1L)
#    checkIdentical(near[[1]], 3L)
#    checkIdentical(near[[2]], 1L)

#    near <- findKNN(g3, g4, ignore.overlap = TRUE)
#    checkIdentical(length(near), 2L)
#    checkIdentical(length(near[[1]]), 0L)
#    checkIdentical(length(near[[2]]), 0L)

    ## ignore.strand, with overlaps in ranges
    near <- findKNN(g1, g2, ignore.strand = TRUE)
    checkIdentical(length(near), 2L)
    checkIdentical(length(near[[1]]), 1L)
    checkIdentical(length(near[[2]]), 1L)
    checkIdentical(near[[1]], 1L)
    checkIdentical(near[[2]], 2L)

    near2 <- findKNN(g3, g4, ignore.strand = TRUE)
    checkIdentical(length(near), length(near2))
    checkIdentical(length(near)[[1]], length(near2)[[1]])
    ## NOTE: I don't know what this test means
#    checkIdentical(length(near)[[2]], length(near2)[[2]])
    checkIdentical(near[[1]], near2[[1]])
    checkIdentical(near[[2]], near2[[2]])

    ## select == "all"
    ## No ties
    checkIdentical(findKNN(g1, g2), findKNN(g1, g2, select = "all"))
    ## With ties
    r <- IRanges(c(14, 1, 6, 10, 10, 14), c(16, 4, 8, 12, 12, 16))
    g <- GRanges("chr1", r, "+")
    target <- findKNN(g, select="all")
    checkIdentical(lengths(target), c(1, 1, 3, 1, 1, 1))
    checkIdentical(target[[3]], c(4, 5, 2))

    ## Case from issue #76
    broads <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    broads <- resize(broads, width(broads)+3000, fix="end")
    gr <- GRanges(
        seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges = IRanges(101:110, end = 111:120, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        score = 1:10,
        GC = seq(1, 0, length=10))
    near <- nearest(gr,broads[ seqnames(broads) %in% seqlevels(gr),])
    target <- findKNN(gr,broads[ seqnames(broads) %in% seqlevels(gr),])
    target10 <- findKNN(gr,broads[ seqnames(broads) %in% seqlevels(gr),], k=10)
    checkIdentical(near, unlist(target))
    checkTrue(all(lengths(target) == 1))
    checkTrue(all(lengths(target10) == 10))
}
