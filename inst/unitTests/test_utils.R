###

test_mergeNamedAtomicVectors <- function()
{
    ## The function is currently not exported.
    mergeNamedAtomicVectors <- GenomicRanges:::mergeNamedAtomicVectors

    x <- c(a=1, b=2, c=3, d=NA)
    checkIdentical(mergeNamedAtomicVectors(x, x), x)

    y1 <- x[FALSE]
    checkIdentical(mergeNamedAtomicVectors(x, y1), x)
    checkIdentical(mergeNamedAtomicVectors(y1, x), x)

    y2 <- c(c=NA_real_, d=NA_real_, b=NA_real_)
    checkIdentical(mergeNamedAtomicVectors(x, y2), x)
    checkIdentical(mergeNamedAtomicVectors(y2, x)[names(x)], x)
    checkIdentical(mergeNamedAtomicVectors(y2, x)[names(y2)], x[names(y2)])

    y3 <- c(c=3, e=NA, d=4, b=2)
    got <- mergeNamedAtomicVectors(x, y3)
    want <- c(a=1, b=2, c=3, d=4, e=NA)
    checkIdentical(got, want)
    got <- mergeNamedAtomicVectors(y3, x)
    want <- c(c=3, e=NA, d=4, b=2, a=1)
    checkIdentical(got, want)

    y4 <- c(c=0, e=NA, d=4, b=0)
    checkException(mergeNamedAtomicVectors(x, y4), silent=TRUE)
    checkException(mergeNamedAtomicVectors(y4, x), silent=TRUE)

    x2 <- c(a=1, b=NA, c=3, d=NA)
    y5 <- c(c=NA, e=5, d=4, b=2)
    got <- mergeNamedAtomicVectors(x2, y5)
    want <- c(a=1, b=2, c=3, d=4, e=5)
    checkIdentical(got, want)
}

