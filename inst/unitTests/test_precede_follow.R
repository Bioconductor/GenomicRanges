test_precede_follow_multiple_ranges <- function()
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
}

test_precede_follow_strand <- function()
{
    subject <- GRanges("chr1", IRanges(rep(1:10, 3), width=1),
                       strand=Rle(strand(c("+", "*", "-")), c(10, 10, 10)))

    ## 'x' on + strand

    x <- GRanges("chr1:4-7:+")

    current <- precede(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(8L, 18L, 28L), sort(subjectHits(current)))
    target <- precede(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "first"), target)
    current <- follow(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(3L, 13L, 23L), sort(subjectHits(current)))
    target <- follow(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "last"), target)

    current <- precede(x, subject, select="all")
    checkIdentical(c(8L, 18L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, subject))
    current <- follow(x, subject, select="all")
    checkIdentical(c(3L, 13L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, subject))

    current <- precede(x, rev(subject), select="all")
    checkIdentical(c(13L, 23L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, rev(subject)))
    current <- follow(x, rev(subject), select="all")
    checkIdentical(c(18L, 28L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, rev(subject)))

    x <- GRanges("chr1:15-20:+")

    current <- precede(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(0L, length(current))
    target <- precede(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "first"), target)
    current <- follow(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(10L, 20L, 30L), sort(subjectHits(current)))
    target <- follow(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "last"), target)

    current <- precede(x, subject, select="all")
    checkIdentical(0L, length(current))
    checkIdentical(selectHits(current, "first"), precede(x, subject))
    current <- follow(x, subject, select="all")
    checkIdentical(c(10L, 20L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, subject))

    current <- precede(x, rev(subject), select="all")
    checkIdentical(0L, length(current))
    checkIdentical(selectHits(current, "first"), precede(x, rev(subject)))
    current <- follow(x, rev(subject), select="all")
    checkIdentical(c(11L, 21L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, rev(subject)))

    ## 'x' on - strand

    x <- GRanges("chr1:4-7:-")

    current <- precede(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(8L, 18L, 28L), sort(subjectHits(current)))
    target <- precede(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "first"), target)
    current <- follow(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(3L, 13L, 23L), sort(subjectHits(current)))
    target <- follow(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "last"), target)

    current <- precede(x, subject, select="all")
    checkIdentical(c(13L, 23L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, subject))
    current <- follow(x, subject, select="all")
    checkIdentical(c(18L, 28L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, subject))

    current <- precede(x, rev(subject), select="all")
    checkIdentical(c(8L, 18L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, rev(subject)))
    current <- follow(x, rev(subject), select="all")
    checkIdentical(c(3L, 13L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, rev(subject)))

    x <- GRanges("chr1:15-20:-")

    current <- precede(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(0L, length(current))
    target <- precede(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "first"), target)
    current <- follow(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(10L, 20L, 30L), sort(subjectHits(current)))
    target <- follow(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "last"), target)

    current <- precede(x, subject, select="all")
    checkIdentical(c(20L, 30L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, subject))
    current <- follow(x, subject, select="all")
    checkIdentical(0L, length(current))
    checkIdentical(selectHits(current, "last"), follow(x, subject))

    current <- precede(x, rev(subject), select="all")
    checkIdentical(c(1L, 11L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, rev(subject)))
    current <- follow(x, rev(subject), select="all")
    checkIdentical(0L, length(current))
    checkIdentical(selectHits(current, "last"), follow(x, rev(subject)))

    ## 'x' on * strand (i.e. on *both* strands)

    x <- GRanges("chr1:4-7")

    current <- precede(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(8L, 18L, 28L), sort(subjectHits(current)))
    target <- precede(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "first"), target)
    current <- follow(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(3L, 13L, 23L), sort(subjectHits(current)))
    target <- follow(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "last"), target)

    current <- precede(x, subject, select="all")
    checkIdentical(c(8L, 18L, 23L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, subject))
    current <- follow(x, subject, select="all")
    checkIdentical(c(3L, 13L, 28L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, subject))

    current <- precede(x, rev(subject), select="all")
    checkIdentical(c(8L, 13L, 23L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, rev(subject)))
    current <- follow(x, rev(subject), select="all")
    checkIdentical(c(3L, 18L, 28L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, rev(subject)))

    x <- GRanges("chr1:15-20")

    current <- precede(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(0L, length(current))
    target <- precede(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "first"), target)
    current <- follow(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(10L, 20L, 30L), sort(subjectHits(current)))
    target <- follow(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "last"), target)

    current <- precede(x, subject, select="all")
    checkIdentical(30L, subjectHits(current))
    checkIdentical(selectHits(current, "first"), precede(x, subject))
    current <- follow(x, subject, select="all")
    checkIdentical(c(10L, 20L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, subject))

    current <- precede(x, rev(subject), select="all")
    checkIdentical(1L, subjectHits(current))
    checkIdentical(selectHits(current, "first"), precede(x, rev(subject)))
    current <- follow(x, rev(subject), select="all")
    checkIdentical(c(11L, 21L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, rev(subject)))

}

test_precede_follow_zero_width_range <- function()
{
    subject <- GRanges("chr1", IRanges(rep(1:10, 3), width=1),
                       strand=Rle(strand(c("+", "*", "-")), c(10, 10, 10)))

    x <- GRanges("chr1:1-0:-")  # zero-width range

    current <- precede(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(1L, 11L, 21L), sort(subjectHits(current)))
    target <- precede(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "first"), target)
    current <- follow(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(0L, length(current))
    target <- follow(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "last"), target)

    current <- precede(x, subject, select="all")
    checkIdentical(0L, length(current))
    checkIdentical(selectHits(current, "first"), precede(x, subject))
    current <- follow(x, subject, select="all")
    checkIdentical(c(11L, 21L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, subject))

    current <- precede(x, rev(subject), select="all")
    checkIdentical(0L, length(current))
    checkIdentical(selectHits(current, "first"), precede(x, rev(subject)))
    current <- follow(x, rev(subject), select="all")
    checkIdentical(c(10L, 20L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "last"), follow(x, rev(subject)))

    x <- GRanges("chr1:1-0:+")  # zero-width range

    current <- precede(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(1L, 11L, 21L), sort(subjectHits(current)))
    target <- precede(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "first"), target)
    current <- follow(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(0L, length(current))
    target <- follow(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "last"), target)

    current <- precede(x, subject, select="all")
    checkIdentical(c(1L, 11L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, subject))
    current <- follow(x, subject, select="all")
    checkIdentical(0L, length(current))
    checkIdentical(selectHits(current, "last"), follow(x, subject))

    current <- precede(x, rev(subject), select="all")
    checkIdentical(c(20L, 30L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, rev(subject)))
    current <- follow(x, rev(subject), select="all")
    checkIdentical(0L, length(current))
    checkIdentical(selectHits(current, "last"), follow(x, rev(subject)))

    x <- GRanges("chr1:1-0")  # zero-width range

    current <- precede(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(c(1L, 11L, 21L), sort(subjectHits(current)))
    target <- precede(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "first"), target)
    current <- follow(x, subject, select="all", ignore.strand=TRUE)
    checkIdentical(0L, length(current))
    target <- follow(x, subject, ignore.strand=TRUE)
    checkIdentical(selectHits(current, "last"), target)

    current <- precede(x, subject, select="all")
    checkIdentical(c(1L, 11L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, subject))
    current <- follow(x, subject, select="all")
    checkIdentical(21L, subjectHits(current))
    checkIdentical(selectHits(current, "last"), follow(x, subject))

    current <- precede(x, rev(subject), select="all")
    checkIdentical(c(20L, 30L), sort(subjectHits(current)))
    checkIdentical(selectHits(current, "first"), precede(x, rev(subject)))
    current <- follow(x, rev(subject), select="all")
    checkIdentical(10L, subjectHits(current))
    checkIdentical(selectHits(current, "last"), follow(x, rev(subject)))
}
