make_subject <- function()
    GRanges(Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
            IRanges(1:10, end=10),
            Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
            seqinfo=Seqinfo(paste0("chr", 1:3)),
            score=1:10, GC=seq(1, 0, length=10))

make_query <- function() {
    GRangesList(
        nomatch =GRanges("chr1", IRanges(start=5, end=10), "+"),
        onematch=GRanges("chr3", IRanges(start=2, end=7), "-"),
        twomatch=GRanges("chr1", IRanges(start=1, end=5), "-"))
}

.checkHits <- function(q_hits, s_hits, q_len, s_len, current, select)
{
    target <- Hits(q_hits, s_hits, q_len, s_len, sort.by.query=TRUE)
    checkIdentical(t(selectHits(target, select=select)), t(unname(current)))
}

test_findOverlaps_no_overlaps_returns_empty_matches <- function()
{
    query <- make_query()
    subject <- make_subject()
    ranges(subject) <- shift(ranges(subject), 1000L)

    ## select = "all"
    for (type in c("any", "start", "end")) {
        current <- findOverlaps(query, subject, type = type, select = "all")
        .checkHits(integer(0), integer(0), 3, 10, current, select="all")

        ans <- countOverlaps(query, subject, type = type)
        checkIdentical(structure(c(0L, 0L, 0L),
                                 names=c("nomatch", "onematch", "twomatch")),
                                 ans)

        ans <- subsetByOverlaps(query, subject, type = type)
        checkIdentical(query[integer(0)], ans)
    }

    ## select = "first"
    expect <- rep(NA_integer_, length(query))
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, type = type, select = "first")
        checkIdentical(expect, ans)
    }
}

test_findOverlaps_empty_query <- function()
{
    query <- GRangesList()
    subject <- make_subject()

    ## select = "all"
    for (type in c("any", "start", "end")) {
        current <- findOverlaps(query, subject, type = type, select = "all")
        .checkHits(integer(0), integer(0), 0, 10, current, select="all")

        ans <- countOverlaps(query, subject, type = type)
        checkIdentical(integer(0), ans)

        ans <- subsetByOverlaps(query, subject, type = type)
        checkIdentical(query, ans)
    }

    ## select = "first"
    expect <- integer()
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, type = type, select = "first")
        checkIdentical(expect, ans)
    }
}

test_findOverlaps_empty_subject <- function()
{
    query <- make_query()
    subject <- GRanges()

    ## select = "all"
    for (type in c("any", "start", "end")) {
        current <- findOverlaps(query, subject, type = type, select = "all")
        .checkHits(integer(0), integer(0), 3, 0, current, select="all")

        ans <- countOverlaps(query, subject, type = type)
        checkIdentical(structure(c(0L, 0L, 0L),
                                 names=c("nomatch", "onematch", "twomatch")),
                       ans)

        ans <- subsetByOverlaps(query, subject, type = type)
        checkIdentical(query[integer(0)], ans)
    }

    ## select = "first"
    expect <- rep(NA_integer_, length(query))
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, type = type, select = "first")
        checkIdentical(expect, ans)
    }
}

test_findOverlaps_zero_one_two_matches <- function()
{
    query <- make_query()
    subject <- make_subject()

    ## select = "all"
    ansAny <- findOverlaps(query, subject, type="any", select="all")
    ansStart <- findOverlaps(query, subject, type="start", select="all")
    ansEnd <- findOverlaps(query, subject, type="end", select="all")
    .checkHits(c(2, 3, 3), c(7, 1, 5), 3, 10, ansAny, select="all")
    .checkHits(3, 1, 3, 10, ansStart, select="all")
    .checkHits(integer(0), integer(0), 3, 10, ansEnd, select="all")

    countsAny <- countOverlaps(query, subject, type="any")
    countsStart <- countOverlaps(query, subject, type="start")
    countsEnd <- countOverlaps(query, subject, type="end")
    .checkHits(c(2, 3, 3), c(7, 1, 5), 3, 10, countsAny, select="count")
    .checkHits(3, 1, 3, 10, countsStart, select="count")
    .checkHits(integer(0), integer(0), 3, 10, countsEnd, select="count")

    subsetAny <- subsetByOverlaps(query, subject, type="any")
    subsetStart <- subsetByOverlaps(query, subject, type="start")
    subsetEnd <- subsetByOverlaps(query, subject, type="end")
    checkIdentical(query[countsAny > 0], subsetAny)
    checkIdentical(query[countsStart > 0], subsetStart)
    checkIdentical(query[countsEnd > 0], subsetEnd)

    ## select = "first"
    ansAny <- findOverlaps(query, subject, type="any", select="first")
    ansStart <- findOverlaps(query, subject, type="start", select="first")
    ansEnd <- findOverlaps(query, subject, type="end", select="first")
    .checkHits(c(2, 3, 3), c(7, 1, 5), 3, 10, ansAny, select="first")
    .checkHits(3, 1, 3, 10, ansStart, select="first")
    .checkHits(integer(0), integer(0), 3, 10, ansEnd, select="first")
}

test_findOverlaps_multimatch_within_one_query <- function()
{
    query <- make_query()
    query[[3L]] <- c(query[[3L]], query[[3L]])
    subject <- make_subject()

    ## select = "all"
    ansAny <- findOverlaps(query, subject, type="any", select="all")
    ansStart <- findOverlaps(query, subject, type="start", select="all")
    ansEnd <- findOverlaps(query, subject, type="end", select="all")
    .checkHits(c(2, 3, 3), c(7, 1, 5), 3, 10, ansAny, select="all")
    .checkHits(3, 1, 3, 10, ansStart, select="all")
    .checkHits(integer(0), integer(0), 3, 10, ansEnd, select="all")

    countsAny <- countOverlaps(query, subject, type="any")
    countsStart <- countOverlaps(query, subject, type="start")
    countsEnd <- countOverlaps(query, subject, type="end")
    .checkHits(c(2, 3, 3), c(7, 1, 5), 3, 10, countsAny, select="count")
    .checkHits(3, 1, 3, 10, countsStart, select="count")
    .checkHits(integer(0), integer(0), 3, 10, countsEnd, select="count")

    subsetAny <- subsetByOverlaps(query, subject, type="any")
    subsetStart <- subsetByOverlaps(query, subject, type="start")
    subsetEnd <- subsetByOverlaps(query, subject, type="end")
    checkIdentical(query[countsAny > 0], subsetAny)
    checkIdentical(query[countsStart > 0], subsetStart)
    checkIdentical(query[countsEnd > 0], subsetEnd)

    ## select = "first"
    ansAny <- findOverlaps(query, subject, type="any", select="first")
    ansStart <- findOverlaps(query, subject, type="start", select="first")
    ansEnd <- findOverlaps(query, subject, type="end", select="first")
    .checkHits(c(2, 3, 3), c(7, 1, 5), 3, 10, ansAny, select="first")
    .checkHits(3, 1, 3, 10, ansStart, select="first")
    .checkHits(integer(0), integer(0), 3, 10, ansEnd, select="first")
}

test_findOverlaps_either_strand <- function()
{
    query <- make_query()
    subject <- make_subject()

    query@unlistData@strand <- Rle(strand(c("*", "*", "-")))

    ## select = "all"
    ansAny <- findOverlaps(query, subject, type="any", select="all")
    ansStart <- findOverlaps(query, subject, type="start", select="all")
    ansEnd <- findOverlaps(query, subject, type="end", select="all")
    .checkHits(c(1, 1, 1, 2, 3, 3), c(1, 5, 6, 7, 1, 5), 3, 10,
               ansAny, select="all")
    .checkHits(c(1, 3), c(5, 1), 3, 10, ansStart, select="all")
    .checkHits(c(1, 1, 1), c(1, 5, 6), 3, 10, ansEnd, select="all")

    countsAny <- countOverlaps(query, subject, type = "any")
    countsStart <- countOverlaps(query, subject, type = "start")
    countsEnd <- countOverlaps(query, subject, type = "end")
    .checkHits(c(1, 1, 1, 2, 3, 3), c(1, 5, 6, 7, 1, 5), 3, 10,
               countsAny, select="count")
    .checkHits(c(1, 3), c(5, 1), 3, 10, countsStart, select="count")
    .checkHits(c(1, 1, 1), c(1, 5, 6), 3, 10, countsEnd, select="count")

    subsetAny <- subsetByOverlaps(query, subject, type = "any")
    subsetStart <- subsetByOverlaps(query, subject, type = "start")
    subsetEnd <- subsetByOverlaps(query, subject, type = "end")
    checkIdentical(query[countsAny > 0], subsetAny)
    checkIdentical(query[countsStart > 0], subsetStart)
    checkIdentical(query[countsEnd > 0], subsetEnd)

    # select = "first"
    ansAny <- findOverlaps(query, subject, type="any", select="first")
    ansStart <- findOverlaps(query, subject, type="start", select="first")
    ansEnd <- findOverlaps(query, subject, type="end", select="first")
    .checkHits(c(1, 1, 1, 2, 3, 3), c(1, 5, 6, 7, 1, 5), 3, 10,
               ansAny, select="first")
    .checkHits(c(1, 3), c(5, 1), 3, 10, ansStart, select="first")
    .checkHits(c(1, 1, 1), c(1, 5, 6), 3, 10, ansEnd, select="first")
}

test_findOverlaps_minoverlap_GRanges_GRangesList <- function()
{
     query <- make_subject()
     subject <- make_query()
     current <- findOverlaps(query, subject, minoverlap = 5)
     .checkHits(1, 3, 10, 3, current, select="all")

     current <- findOverlaps(query, subject, minoverlap = 6)
     .checkHits(integer(0), integer(0), 10, 3, current, select="all")
}

test_findOverlaps_minoverlap_GRangesList_GRanges <- function()
{
     subject <- make_subject()
     query <- make_query()
     current <- findOverlaps(query, subject, minoverlap = 5)
     .checkHits(3, 1, 3, 10, current, select="all")

     current <- findOverlaps(query, subject, minoverlap = 6)
     .checkHits(integer(0), integer(0), 3, 10, current, select="all")
}

test_findOverlaps_minoverlap_GRangesList_GRangesList <- function()
{
     query <- make_query()
     subject <- GRangesList("g1" = make_subject())
     current <- findOverlaps(query, subject, minoverlap = 1)
     .checkHits(c(2, 3), c(1, 1), 3, 1, current, select="all")

     query <- make_query()
     subject <- GRangesList("g1" = make_subject())
     current <- findOverlaps(query, subject, minoverlap = 6)
     .checkHits(3, 1, 3, 1, current, select="all")

     query <- make_query()
     subject <- GRangesList("g1" = make_subject())
     current <- findOverlaps(query, subject, minoverlap = 7)
     .checkHits(integer(0), integer(0), 3, 1, current, select="all")

     current <- findOverlaps(subject, query, minoverlap = 6)
     .checkHits(1, 3, 1, 3, current, select="all")
}

test_findOverlaps_with_circular_sequences <- function()
{
    gr <- GRanges(seqnames=rep.int("A", 4),
                  ranges=IRanges(start=c(2, 4, 6, 8), width=3))

    ## With A of length 9 --> no overlap between last and first ranges.
    gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=9, isCircular=TRUE)
    current0 <- findOverlaps(gr, gr)
    target0_q_hits <- c(1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L)
    target0_s_hits <- c(1L, 2L, 1L, 2L, 3L, 2L, 3L, 4L, 3L, 4L)
    .checkHits(target0_q_hits, target0_s_hits, 4, 4, current0, select="all")

    ## With A of length 8 --> last and first ranges do overlap.
    gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=8, isCircular=TRUE)
    current1 <- findOverlaps(gr, gr)
    .checkHits(c(1, target0_q_hits, 4), c(4, target0_s_hits, 1), 4, 4,
               current1, select="all")

    ## With A of length 8 and minoverlap=2 --> no overlap between last
    ## and first ranges.
    current2 <- findOverlaps(gr, gr, minoverlap=2)
    .checkHits(1:4, 1:4, 4, 4, current2, select="all")

    ## With A of length 7 and minoverlap=2 --> last and first ranges
    ## do overlap.
    gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=7, isCircular=TRUE)
    current3 <- findOverlaps(gr, gr, minoverlap=2)
    .checkHits(c(1, 1:4, 4), c(4, 1:4, 1), 4, 4, current3, select="all")

    ## type = "within"
    q0 <- GRanges("A", IRanges(c(11, 5, 4, 11, 11, 4),
                     c(30, 30, 30, 50, 51, 51)))
    s0 <- GRanges("A", IRanges(5, width=46))
    s0@seqinfo <- Seqinfo(seqnames="A", seqlengths=100, isCircular=TRUE)
    ## sanity check with linear shift
    fo0 <- findOverlaps(q0, s0, type="within")
    expected <- c(1L, 2L, 4L)
    checkIdentical(queryHits(fo0), expected)
    A=90
    q1 <- shift(q0, A)
    s1 <- shift(s0, A)
    fo1 <- findOverlaps(q1, s1, type="within")
    checkIdentical(queryHits(fo1), expected)
    ## circular shift
    n1=-1; n2=0
    q2 <- shift(q0, A + 100 * n1)
    s2 <- shift(s0, A + 100 * n2)
    fo1 <- findOverlaps(q1, s1, type="within")
    checkIdentical(queryHits(fo1), expected)

    ## With A of length 8 --> range 3 is within range 2
    gr <- GRanges(seqnames=rep.int("A", 4),
                  ranges=IRanges(start=c(2, 4, 6, 8), width=c(3, 3, 3, 5)))
    gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=8, isCircular=TRUE)
    current4 <- findOverlaps(gr, gr, type="within")
    .checkHits(c(1, 1:4), c(1, 4, 2, 3, 4), 4, 4, current4, select="all")

    ## With A of length 9 --> range 3 is not within range 2
    gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=9, isCircular=TRUE)
    current5 <- findOverlaps(gr, gr, type="within")
    .checkHits(1:4, 1:4, 4, 4, current5, select="all")
}

test_poverlaps <- function() {
    ans <- poverlaps(GRanges(), GRanges())
    checkIdentical(ans, Rle())

    ans <- poverlaps(GRanges("chr1:11-15"), GRanges("chr1:16-20"))
    checkIdentical(ans, Rle(FALSE))

    ans <- poverlaps(GRanges("chr1:11-16"), GRanges("chr1:16-20"))
    checkIdentical(ans, Rle(TRUE))

    ans <- poverlaps(GRanges(c("chr1:11-15", "chr1:11-16")),
                     GRanges("chr1:16-20"))
    checkIdentical(ans, Rle(c(FALSE, TRUE)))
}

test_findOverlaps_character <- function() {
    gr <- GRanges(c("chr1:11-15", "chr1:5-100", "chr2:10-11"))

    # GRanges as the query, character as the subject.
    {
        subj <- "chr1"
        ans <- findOverlaps(gr, subj)
        checkIdentical(queryHits(ans), c(1L, 2L))
        checkIdentical(subjectHits(ans), c(1L, 1L))
        checkIdentical(countOverlaps(gr, subj), c(1L, 1L, 0L))
        checkIdentical(overlapsAny(gr, subj), c(TRUE, TRUE, FALSE))
        checkIdentical(ans, findOverlaps(gr, subj, type="within"))
        checkIdentical(subsetByOverlaps(gr, subj), gr[1:2])
        checkIdentical(subsetByOverlaps(gr, subj, invert=TRUE), gr[3])

        # Trying with duplicates and missing chromosomes.
        subj <- c("chr2", "chr1", "chr3", "chr2")
        ans <- findOverlaps(gr, subj)
        checkIdentical(queryHits(ans), c(1L, 2L, 3L, 3L))
        checkIdentical(subjectHits(ans), c(2L, 2L, 1L, 4L))
        checkIdentical(countOverlaps(gr, subj), c(1L, 1L, 2L))
        checkIdentical(overlapsAny(gr, subj), c(TRUE, TRUE, TRUE))
        checkIdentical(ans, findOverlaps(gr, subj, type="within"))
        checkIdentical(subsetByOverlaps(gr, subj), gr)
        checkIdentical(subsetByOverlaps(gr, subj, invert=TRUE), gr[integer()])

        checkException(findOverlaps(gr, "chr2", type="start"), silent=TRUE)
        checkException(findOverlaps(gr, "chr2", minoverlap=1), silent=TRUE)
    }

    # character as the query, GRanges as the subject.
    {
        qry <- "chr2"
        ans <- findOverlaps(qry, gr)
        checkIdentical(queryHits(ans), 1L)
        checkIdentical(subjectHits(ans), 3L)
        checkIdentical(countOverlaps(qry, gr), 1L)
        checkIdentical(overlapsAny(qry, gr), TRUE)
        checkIdentical(subsetByOverlaps(qry, gr), qry)
        checkIdentical(subsetByOverlaps(qry, gr, invert=TRUE), character(0))

        # Trying with duplicates and missing chromosomes.
        qry <- c("chr3", "chr1", "chr2", "chr1")
        ans <- findOverlaps(qry, gr)
        checkIdentical(queryHits(ans), c(2L, 2L, 3L, 4L, 4L))
        checkIdentical(subjectHits(ans), c(1L, 2L, 3L, 1L, 2L))
        checkIdentical(countOverlaps(qry, gr), c(0L, 2L, 1L, 2L))
        checkIdentical(overlapsAny(qry, gr), c(FALSE, TRUE, TRUE, TRUE))
        checkIdentical(subsetByOverlaps(qry, gr), qry[2:4])
        checkIdentical(subsetByOverlaps(qry, gr, invert=TRUE), qry[1])

        checkException(findOverlaps("chr2", gr, type="within"), silent=TRUE)
        checkException(findOverlaps("chr2", gr, minoverlap=1), silent=TRUE) 
    }

    gr <- GRanges(c("chr1:11-15", "chr1:5-100", "chr2:1-10", "chr2:10-11", "chr2:5-7"))
    grl <- unname(split(gr, c(1,2,2,3,3)))

    # GRangesList as the query, character as the subject.
    {
        subj <- "chr1"
        ans <- findOverlaps(grl, subj)
        checkIdentical(queryHits(ans), c(1L, 2L))
        checkIdentical(subjectHits(ans), c(1L, 1L))
        checkIdentical(countOverlaps(grl, subj), c(1L, 1L, 0L))
        checkIdentical(overlapsAny(grl, subj), c(TRUE, TRUE, FALSE))
        checkIdentical(subsetByOverlaps(grl, subj), grl[1:2])
        checkIdentical(subsetByOverlaps(grl, subj, invert=TRUE), grl[3])

        ans <- findOverlaps(grl, subj, type="within")
        checkIdentical(queryHits(ans), 1L) 
        checkIdentical(subjectHits(ans), 1L)
        checkIdentical(countOverlaps(grl, subj, type="within"), c(1L, 0L, 0L))
        checkIdentical(overlapsAny(grl, subj, type="within"), c(TRUE, FALSE, FALSE))
        checkIdentical(subsetByOverlaps(grl, subj, type="within"), grl[1])
        checkIdentical(subsetByOverlaps(grl, subj, type="within", invert=TRUE), grl[2:3])

        # Trying with duplicates and missing chromosomes.
        subj <- c("chr3", "chr2", "chr1", "chr1")
        ans <- sort(findOverlaps(grl, subj)) # sorting for an easier comparison to the expected result.
        checkIdentical(queryHits(ans), c(1L, 1L, 2L, 2L, 2L, 3L))
        checkIdentical(subjectHits(ans), c(3L, 4L, 2L, 3L, 4L, 2L))
        checkIdentical(countOverlaps(grl, subj), c(2L, 3L, 1L))
        checkIdentical(overlapsAny(grl, subj), c(TRUE, TRUE, TRUE))
        checkIdentical(subsetByOverlaps(grl, subj), grl)
        checkIdentical(subsetByOverlaps(grl, subj, invert=TRUE), grl[integer(0)])

        ans <- findOverlaps(grl, subj, type="within")
        checkIdentical(queryHits(ans), c(1L, 1L, 3L)) 
        checkIdentical(subjectHits(ans), c(3L, 4L, 2L))
        checkIdentical(countOverlaps(grl, subj, type="within"), c(2L, 0L, 1L))
        checkIdentical(overlapsAny(grl, subj, type="within"), c(TRUE, FALSE, TRUE))
        checkIdentical(subsetByOverlaps(grl, subj, type="within"), grl[c(1,3)])
        checkIdentical(subsetByOverlaps(grl, subj, type="within", invert=TRUE), grl[2])

        checkException(findOverlaps(grl, "chr2", type="start"), silent=TRUE)
        checkException(findOverlaps(grl, "chr2", minoverlap=1), silent=TRUE)
    }

    # character as the query, GRangesList as the subject.
    {
        qry <- "chr2"
        ans <- findOverlaps(qry, grl)
        checkIdentical(queryHits(ans), c(1L, 1L))
        checkIdentical(subjectHits(ans), c(2L, 3L))
        checkIdentical(countOverlaps(qry, grl), 2L)
        checkIdentical(overlapsAny(qry, grl), TRUE)
        checkIdentical(subsetByOverlaps(qry, grl), qry)
        checkIdentical(subsetByOverlaps(qry, grl, invert=TRUE), character(0))

        # Trying with duplicates and missing chromosomes.
        qry <- c("chr1", "chr3", "chr2", "chr1")
        ans <- findOverlaps(qry, grl)
        checkIdentical(queryHits(ans), c(1L, 1L, 3L, 3L, 4L, 4L))
        checkIdentical(subjectHits(ans), c(1L, 2L, 2L, 3L, 1L, 2L))
        checkIdentical(countOverlaps(qry, grl), c(2L, 0L, 2L, 2L))
        checkIdentical(overlapsAny(qry, grl), c(TRUE, FALSE, TRUE, TRUE))
        checkIdentical(subsetByOverlaps(qry, grl), qry[c(1L, 3L, 4L)])
        checkIdentical(subsetByOverlaps(qry, grl, invert=TRUE), qry[2])

        checkException(findOverlaps("chr2", grl, type="within"), silent=TRUE)
        checkException(findOverlaps("chr2", grl, minoverlap=1), silent=TRUE)
    }
}
