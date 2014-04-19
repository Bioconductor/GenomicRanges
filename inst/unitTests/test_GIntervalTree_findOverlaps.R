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
    subject <- GIntervalTree(subject)

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
    subject <- GIntervalTree(make_subject())

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
    subject <- GIntervalTree(new("GRanges"))

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
    subject <- GIntervalTree(make_subject())

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
    subject <- GIntervalTree(make_subject())

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
    subject <- GIntervalTree(make_subject())

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

# this needs GIntervalTreeList?
# test_findOverlaps_minoverlap_GRanges_GRangesList <- function() {
    
#     DEACTIVATED()
#      query <- make_subject()
#      subject <- make_query()
#      current <- findOverlaps(query, subject, minoverlap = 5)
#      target <-  new("Hits",
#                     queryHits = 1L, subjectHits = 3L,
#                     queryLength = 10L, subjectLength = 3L)
#      checkIdentical(target, current)

#      current <- findOverlaps(query, subject, minoverlap = 6)
#      target <-  new("Hits",
#                     queryHits = integer(0), subjectHits = integer(0),
#                     queryLength = 10L, subjectLength = 3L)
#      checkIdentical(target, current)
# }


test_findOverlaps_minoverlap_GRangesList_GRanges <- function() {
    
     subject <- GIntervalTree(make_subject())
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

test_findOverlaps_seqlevelIssues <- function() {
    gr=GRanges(seqnames=c("chr1","chr2"),ranges=IRanges(1:2,width=1))
    gtree=GIntervalTree(gr)

    olaps1=findOverlaps(gr,gr[1])
    olaps2=findOverlaps(gr,gtree[1])
    checkIdentical(olaps1, olaps2)

    gr1 = GRanges(seqnames = c("chr1", "chr2"), ranges = IRanges(1:2, width = 1))    
    gr2 = gr1
    seqlevels(gr2, force = TRUE) = c("chr2", "chr1")

    olaps1=findOverlaps(gr1,gr2)
    olaps2=findOverlaps(gr1,GIntervalTree(gr2))
    checkIdentical(olaps1,olaps2)


}
# this needs GIntervalTreeList?
# test_findOverlaps_minoverlap_GrangesList_GRangesList <- function() {
#      query <- make_query()
#      subject <- GRangesList("g1" = make_subject())
#      current <- findOverlaps(query, subject, minoverlap = 1)
#      target <- new("Hits",
#                    queryHits = c(2L, 3L), subjectHits = c(1L, 1L),
#                    queryLength = 3L, subjectLength = 1L)
#      checkIdentical(target, current)
     
#      query <- make_query()
#      subject <- GRangesList("g1" = make_subject())
#      current <- findOverlaps(query, subject, minoverlap = 6)
#      target <- new("Hits",
#                    queryHits = 3L, subjectHits = 1L,
#                    queryLength = 3L, subjectLength = 1L)
#      checkIdentical(target, current)

#      query <- make_query()
#      subject <- GRangesList("g1" = make_subject())
#      current <- findOverlaps(query, subject, minoverlap = 7)
#      target <-  new("Hits",
#                     queryHits = integer(0), subjectHits = integer(0),
#                     queryLength = 3L, subjectLength = 1L)
#      checkIdentical(target, current)

#      current <- findOverlaps(subject, query, minoverlap = 6)
#      target <-  new("Hits",
#                     queryHits = 1L, subjectHits = 3L,
#                     queryLength = 1L, subjectLength = 3L)
#      checkIdentical(target, current)

# }

# this should raise an error
test_findOverlaps_with_circular_sequences <- function()
{
    gr <- GIntervalTree(GRanges(seqnames=rep.int("A", 4),
                  ranges=IRanges(start=c(2, 4, 6, 8), width=3)))

    ## With A of length 9 --> no overlap between last and first ranges.
    gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=9, isCircular=TRUE)
    checkException(current0 <- findOverlaps(gr, gr))
    # matchMatrix0 <- matrix(
    #                     c(1L, 1L, 2L, 2L, 2L, 3L, 3L, 3L, 4L, 4L,
    #                       1L, 2L, 1L, 2L, 3L, 2L, 3L, 4L, 3L, 4L),
    #                     ncol = 2,
    #                     dimnames = list(NULL, c("queryHits", "subjectHits")))
    # target0 <-  new("Hits",
    #                 queryHits = unname(matchMatrix0[ , 1L]),
    #                 subjectHits = unname(matchMatrix0[ , 2L]),
    #                 queryLength = 4L, subjectLength = 4L)
    # checkIdentical(target0, current0)

    # ## With A of length 8 --> last and first ranges do overlap.
    # gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=8, isCircular=TRUE)
    # current1 <- findOverlaps(gr, gr)
    # matchMatrix1 <- rbind(matchMatrix0, matrix(c(1L, 4L, 4L, 1L), ncol = 2))
    # o1 <- S4Vectors:::orderIntegerPairs(matchMatrix1[ , 1L],
    #                                     matchMatrix1[ , 2L])
    # matchMatrix1 <- matchMatrix1[o1, ]
    # target1 <- new("Hits",
    #                queryHits = unname(matchMatrix1[ , 1L]),
    #                subjectHits = unname(matchMatrix1[ , 2L]),
    #                queryLength = 4L, subjectLength = 4L)
    # checkIdentical(target1, current1)

    # ## With A of length 8 and minoverlap=2 --> no overlap between last
    # ## and first ranges.
    # current2 <- findOverlaps(gr, gr, minoverlap=2)
    # matchMatrix2 <- matrix(c(1:4, 1:4), ncol = 2,
    #                        dimnames = list(NULL, c("queryHits", "subjectHits")))
    # target2 <- new("Hits",
    #                queryHits = unname(matchMatrix2[ , 1L]),
    #                subjectHits = unname(matchMatrix2[ , 2L]),
    #                queryLength = 4L, subjectLength = 4L)
    # checkIdentical(target2, current2)

    # ## With A of length 7 and minoverlap=2 --> last and first ranges
    # ## do overlap.
    # gr@seqinfo <- Seqinfo(seqnames="A", seqlengths=7, isCircular=TRUE)
    # current3 <- findOverlaps(gr, gr, minoverlap=2)
    # matchMatrix3 <- rbind(matchMatrix2, matrix(c(1L, 4L, 4L, 1L), ncol = 2))
    # o3 <- S4Vectors:::orderIntegerPairs(matchMatrix3[ , 1L],
    #                                     matchMatrix3[ , 2L])
    # matchMatrix3 <- matchMatrix3[o3, ]
    # target3 <- new("Hits",
    #                queryHits = unname(matchMatrix3[ , 1L]),
    #                subjectHits = unname(matchMatrix3[ , 2L]),
    #                queryLength = 4L, subjectLength = 4L)
    # checkIdentical(target3, current3)
}

