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
    x <- GRanges("chr1", IRanges(c(0, 10, 99), width=1))
    subject <- GRanges("chr1", IRanges(15, width=1))
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
