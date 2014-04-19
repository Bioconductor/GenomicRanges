make_subject <- function() {
    new("GRanges",
        seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1),
        strand = Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
        seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")),
        elementMetadata = DataFrame(score = 1:10, GC = seq(1, 0, length=10)))
}

make_query <- function() {
    GRangesList(nomatch = GRanges(seqnames = "chr1",
                                  ranges = IRanges(start=5, end=10),
                                  strand = "+"),
                onematch = GRanges(seqnames = "chr3",
                                   ranges = IRanges(start=2, end=7),
                                   strand = "-"),
                twomatch = GRanges(seqnames = "chr1",
                                   ranges = IRanges(start=1, end=5),
                                   strand = "-"))
}

test_findOverlaps_no_overlaps_returns_empty_matches <- function()
{
    query <- make_query()
    subject <- make_subject()
    ranges(subject) <- shift(ranges(subject), 1000L)

    ## select = "all"
    expect <- new("Hits", queryHits = integer(0), subjectHits = integer(0),
                          queryLength = 3L, subjectLength = 10L)
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, type = type, select = "all")
        checkIdentical(expect, ans)

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
    query <- new("GRangesList")
    subject <- make_subject()

    ## select = "all"
    expect <- new("Hits", queryHits = integer(0), subjectHits = integer(0),
                          queryLength = 0L, subjectLength = 10L)
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, type = type, select = "all")
        checkIdentical(expect, ans)

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
    subject <- new("GRanges")

    ## select = "all"
    expect <- new("Hits", queryHits = integer(0), subjectHits = integer(0),
                          queryLength = 3L, subjectLength = 0L)
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, type = type, select = "all")
        checkIdentical(expect, ans)

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
    expectAny <- new("Hits",
                     queryHits = c(2L, 3L, 3L), subjectHits = c(7L, 1L, 5L),
                     queryLength = 3L, subjectLength = 10L)
    expectStart <- new("Hits",
                       queryHits = 3L, subjectHits = 1L,
                       queryLength = 3L, subjectLength = 10L)
    expectEnd <- new("Hits",
                     queryHits = integer(0), subjectHits = integer(0),
                     queryLength = 3L, subjectLength = 10L)
    ansAny <- findOverlaps(query, subject, select = "all", type = "any")
    ansStart <- findOverlaps(query, subject, select = "all", type = "start")
    ansEnd <- findOverlaps(query, subject, select = "all", type = "end")
    checkIdentical(expectAny, ansAny)
    checkIdentical(expectStart, ansStart)
    checkIdentical(expectEnd, ansEnd)

    countsAny <- countOverlaps(query, subject, type = "any")
    countsStart <- countOverlaps(query, subject, type = "start")
    countsEnd <- countOverlaps(query, subject, type = "end")
    checkIdentical(structure(tabulate(queryHits(expectAny), 3),
                             names=c("nomatch", "onematch", "twomatch")),
                   countsAny)
    checkIdentical(structure(tabulate(queryHits(expectStart), 3),
                             names=c("nomatch", "onematch", "twomatch")),
                   countsStart)
    checkIdentical(structure(tabulate(queryHits(expectEnd), 3),
                             names=c("nomatch", "onematch", "twomatch")),
                   countsEnd)

    subsetAny <- subsetByOverlaps(query, subject, type = "any")
    subsetStart <- subsetByOverlaps(query, subject, type = "start")
    subsetEnd <- subsetByOverlaps(query, subject, type = "end")
    checkIdentical(query[countsAny > 0], subsetAny)
    checkIdentical(query[countsStart > 0], subsetStart)
    checkIdentical(query[countsEnd > 0], subsetEnd)

    ## select = "first"
    expectAny <- c(NA_integer_, 7L, 1L)
    expectStart <- c(NA_integer_, NA_integer_, 1L)
    expectEnd <- c(NA_integer_, NA_integer_, NA_integer_)
    ansAny <- findOverlaps(query, subject, type = "any", select = "first")
    ansStart <- findOverlaps(query, subject, type = "start", select = "first")
    ansEnd <- findOverlaps(query, subject, type = "end", select = "first")
    checkIdentical(expectAny, ansAny)
    checkIdentical(expectAny, findOverlaps(query, subject, select="first"))
    checkIdentical(expectStart, ansStart)
    checkIdentical(expectEnd, ansEnd)
}

test_findOverlaps_multimatch_within_one_query <- function()
{
    query <- make_query()
    query[[3L]] <- c(query[[3L]], query[[3L]])
    subject <- make_subject()

    ## select = "all"
    expectAny <- new("Hits",
                     queryHits = c(2L, 3L, 3L), subjectHits = c(7L, 1L, 5L),
                     queryLength = 3L, subjectLength = 10L)
    expectStart <- new("Hits",
                       queryHits = 3L, subjectHits = 1L,
                       queryLength = 3L, subjectLength = 10L)
    expectEnd <- new("Hits",
                     queryHits = integer(0), subjectHits = integer(0),
                     queryLength = 3L, subjectLength = 10L)
    ansAny <- findOverlaps(query, subject, select = "all", type = "any")
    ansStart <- findOverlaps(query, subject, select = "all", type = "start")
    ansEnd <- findOverlaps(query, subject, select = "all", type = "end")
    checkIdentical(expectAny, ansAny)
    checkIdentical(expectStart, ansStart)
    checkIdentical(expectEnd, ansEnd)

    countsAny <- countOverlaps(query, subject, type = "any")
    countsStart <- countOverlaps(query, subject, type = "start")
    countsEnd <- countOverlaps(query, subject, type = "end")
    checkIdentical(structure(tabulate(queryHits(expectAny), 3),
                             names=c("nomatch", "onematch", "twomatch")),
                   countsAny)
    checkIdentical(structure(tabulate(queryHits(expectStart), 3),
                             names=c("nomatch", "onematch", "twomatch")),
                   countsStart)
    checkIdentical(structure(tabulate(queryHits(expectEnd), 3),
                             names=c("nomatch", "onematch", "twomatch")),
                   countsEnd)

    subsetAny <- subsetByOverlaps(query, subject, type = "any")
    subsetStart <- subsetByOverlaps(query, subject, type = "start")
    subsetEnd <- subsetByOverlaps(query, subject, type = "end")
    checkIdentical(query[countsAny > 0], subsetAny)
    checkIdentical(query[countsStart > 0], subsetStart)
    checkIdentical(query[countsEnd > 0], subsetEnd)

    ## select = "first"
    expectAny <- c(NA_integer_, 7L, 1L)
    expectStart <- c(NA_integer_, NA_integer_, 1L)
    expectEnd <- c(NA_integer_, NA_integer_, NA_integer_)
    ansAny <- findOverlaps(query, subject, type = "any", select = "first")
    ansStart <- findOverlaps(query, subject, type = "start", select = "first")
    ansEnd <- findOverlaps(query, subject, type = "end", select = "first")
    checkIdentical(expectAny, ansAny)
    checkIdentical(expectAny, findOverlaps(query, subject, select="first"))
    checkIdentical(expectStart, ansStart)
    checkIdentical(expectEnd, ansEnd)
}

test_findOverlaps_either_strand <- function()
{
    query <- make_query()
    subject <- make_subject()

    query@unlistData@strand <- Rle(strand(c("*", "*", "-")))

    ## select = "all"
    expectAny <- new("Hits",
                     queryHits = c(1L, 1L, 1L, 2L, 3L, 3L),
                     subjectHits = c(1L, 5L, 6L, 7L, 1L, 5L),
                     queryLength = 3L, subjectLength = 10L)
    expectStart <- new("Hits",
                       queryHits = c(1L, 3L), subjectHits = c(5L, 1L),
                       queryLength = 3L, subjectLength = 10L)
    expectEnd <- new("Hits",
                     queryHits = c(1L, 1L, 1L), subjectHits = c(1L, 5L, 6L),
                     queryLength = 3L, subjectLength = 10L)
    ansAny <- findOverlaps(query, subject, type = "any", select = "all")
    ansStart <- findOverlaps(query, subject, type = "start", select = "all")
    ansEnd <- findOverlaps(query, subject, type = "end", select = "all")
    checkIdentical(expectAny, ansAny)
    checkIdentical(expectStart, ansStart)
    checkIdentical(expectEnd, ansEnd)

    countsAny <- countOverlaps(query, subject, type = "any")
    countsStart <- countOverlaps(query, subject, type = "start")
    countsEnd <- countOverlaps(query, subject, type = "end")
    checkIdentical(structure(tabulate(queryHits(expectAny), 3),
                             names=c("nomatch", "onematch", "twomatch")),
                   countsAny)
    checkIdentical(structure(tabulate(queryHits(expectStart), 3),
                             names=c("nomatch", "onematch", "twomatch")),
                   countsStart)
    checkIdentical(structure(tabulate(queryHits(expectEnd), 3),
                             names=c("nomatch", "onematch", "twomatch")),
                   countsEnd)

    subsetAny <- subsetByOverlaps(query, subject, type = "any")
    subsetStart <- subsetByOverlaps(query, subject, type = "start")
    subsetEnd <- subsetByOverlaps(query, subject, type = "end")
    checkIdentical(query[countsAny > 0], subsetAny)
    checkIdentical(query[countsStart > 0], subsetStart)
    checkIdentical(query[countsEnd > 0], subsetEnd)

    # select = "first"
    expectAny <- c(1L, 7L, 1L)
    expectStart <- c(5L, NA_integer_, 1L)
    expectEnd <- c(1L, NA_integer_, NA_integer_)
    ansAny <- findOverlaps(query, subject, type = "any", select = "first")
    ansStart <- findOverlaps(query, subject, type = "start", select = "first")
    ansEnd <- findOverlaps(query, subject, type = "end", select = "first")
    checkIdentical(expectAny, ansAny)
    checkIdentical(expectAny, findOverlaps(query, subject, select="first"))
    checkIdentical(expectStart, ansStart)
    checkIdentical(expectEnd, ansEnd)
}

test_findOverlaps_minoverlap_GRanges_GRangesList <- function() {
    
     query <- make_subject()
     subject <- make_query()
     current <- findOverlaps(query, subject, minoverlap = 5)
     target <-  new("Hits",
                    queryHits = 1L, subjectHits = 3L,
                    queryLength = 10L, subjectLength = 3L)
     checkIdentical(target, current)

     current <- findOverlaps(query, subject, minoverlap = 6)
     target <-  new("Hits",
                    queryHits = integer(0), subjectHits = integer(0),
                    queryLength = 10L, subjectLength = 3L)
     checkIdentical(target, current)
}


test_findOverlaps_minoverlap_GRangesList_GRanges <- function() {
    
     subject <- make_subject()
     query <- make_query()
     current <- findOverlaps(query, subject, minoverlap = 5)
     target <-  new("Hits",
                    queryHits = 3L, subjectHits = 1L,
                    queryLength = 3L, subjectLength = 10L)
     checkIdentical(target, current)

     current <- findOverlaps(query, subject, minoverlap = 6)
     target <-  new("Hits",
                    queryHits = integer(0), subjectHits = integer(0),
                    queryLength = 3L, subjectLength = 10L)
     checkIdentical(target, current)
}


test_findOverlaps_minoverlap_GrangesList_GRangesList <- function() {

     query <- make_query()
     subject <- GRangesList("g1" = make_subject())
     current <- findOverlaps(query, subject, minoverlap = 1)
     target <- new("Hits",
                   queryHits = c(2L, 3L), subjectHits = c(1L, 1L),
                   queryLength = 3L, subjectLength = 1L)
     checkIdentical(target, current)
     
     query <- make_query()
     subject <- GRangesList("g1" = make_subject())
     current <- findOverlaps(query, subject, minoverlap = 6)
     target <- new("Hits",
                   queryHits = 3L, subjectHits = 1L,
                   queryLength = 3L, subjectLength = 1L)
     checkIdentical(target, current)

     query <- make_query()
     subject <- GRangesList("g1" = make_subject())
     current <- findOverlaps(query, subject, minoverlap = 7)
     target <-  new("Hits",
                    queryHits = integer(0), subjectHits = integer(0),
                    queryLength = 3L, subjectLength = 1L)
     checkIdentical(target, current)

     current <- findOverlaps(subject, query, minoverlap = 6)
     target <-  new("Hits",
                    queryHits = 1L, subjectHits = 3L,
                    queryLength = 1L, subjectLength = 3L)
     checkIdentical(target, current)

}

test_findOverlaps_with_circular_sequences <- function()
{
    gr <- GRanges(seqnames=rep.int("A", 4),
                  ranges=IRanges(start=c(2, 4, 6, 8), width=3))

    ## With A of length 9 --> no overlap between last and first ranges.
    gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=9, isCircular=TRUE)
    current0 <- findOverlaps(gr, gr)
    matchMatrix0 <- matrix(
                        c(1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L,
                          1L, 2L, 1L, 2L, 3L, 2L, 3L, 4L, 3L, 4L),
                        ncol = 2,
                        dimnames = list(NULL, c("queryHits", "subjectHits")))
    target0 <-  new("Hits",
                    queryHits = unname(matchMatrix0[ , 1L]),
                    subjectHits = unname(matchMatrix0[ , 2L]),
                    queryLength = 4L, subjectLength = 4L)
    checkIdentical(target0, current0)

    ## With A of length 8 --> last and first ranges do overlap.
    gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=8, isCircular=TRUE)
    current1 <- findOverlaps(gr, gr)
    matchMatrix1 <- rbind(matchMatrix0, matrix(c(1L, 4L, 4L, 1L), ncol = 2))
    o1 <- S4Vectors:::orderIntegerPairs(matchMatrix1[ , 1L],
                                        matchMatrix1[ , 2L])
    matchMatrix1 <- matchMatrix1[o1, ]
    target1 <- new("Hits",
                   queryHits = unname(matchMatrix1[ , 1L]),
                   subjectHits = unname(matchMatrix1[ , 2L]),
                   queryLength = 4L, subjectLength = 4L)
    checkIdentical(target1, current1)

    ## With A of length 8 and minoverlap=2 --> no overlap between last
    ## and first ranges.
    current2 <- findOverlaps(gr, gr, minoverlap=2)
    matchMatrix2 <- matrix(c(1:4, 1:4), ncol = 2,
                           dimnames = list(NULL, c("queryHits", "subjectHits")))
    target2 <- new("Hits",
                   queryHits = unname(matchMatrix2[ , 1L]),
                   subjectHits = unname(matchMatrix2[ , 2L]),
                   queryLength = 4L, subjectLength = 4L)
    checkIdentical(target2, current2)

    ## With A of length 7 and minoverlap=2 --> last and first ranges
    ## do overlap.
    gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=7, isCircular=TRUE)
    current3 <- findOverlaps(gr, gr, minoverlap=2)
    matchMatrix3 <- rbind(matchMatrix2, matrix(c(1L, 4L, 4L, 1L), ncol = 2))
    o3 <- S4Vectors:::orderIntegerPairs(matchMatrix3[ , 1L],
                                        matchMatrix3[ , 2L])
    matchMatrix3 <- matchMatrix3[o3, ]
    target3 <- new("Hits",
                   queryHits = unname(matchMatrix3[ , 1L]),
                   subjectHits = unname(matchMatrix3[ , 2L]),
                   queryLength = 4L, subjectLength = 4L)
    checkIdentical(target3, current3)

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
    matchMatrix4 <- matrix(c(1L, 1L, 2L, 3L, 4L, 
                             1L, 4L, 2L, 3L, 4L), ncol = 2,
                           dimnames = list(NULL, c("queryHits", "subjectHits")))
    target4 <- new("Hits",
                   queryHits = unname(matchMatrix4[ , 1L]),
                   subjectHits = unname(matchMatrix4[ , 2L]),
                   queryLength = 4L, subjectLength = 4L)
    checkIdentical(target4, current4)
    ## With A of length 9 --> range 3 is not within range 2 
    gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=9, isCircular=TRUE)
    current5 <- findOverlaps(gr, gr, type="within")
    target5 <- new("Hits",
                   queryHits = unname(matchMatrix2[ , 1L]),
                   subjectHits = unname(matchMatrix2[ , 2L]),
                   queryLength = 4L, subjectLength = 4L)
    checkIdentical(target5, current5)
}
