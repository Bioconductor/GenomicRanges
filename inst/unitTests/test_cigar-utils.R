###

test_cigarRangesAlongReferenceSpace <- function()
{
    cigar <- c("30M5000N10M", "50M4S", "90=10X5I50M10D40M", "18M10I22M", "99I")
    pos <- c(101, 201, 1001,  301, 2001)
    seqlevels <- c("chr2", "chr6")
    rname <- factor(c("chr6", "chr6", "chr2", "chr6", "chr2"),
                    levels=seqlevels)

    ops <- c("M", "=", "X", "I", "D")
    current <- as.list(cigarRangesAlongReferenceSpace(cigar, pos=pos, f=rname,
                                                      ops=ops))
    ir2 <- IRanges(start=c(1001, 1091, 1101, 1101, 1151, 1161, 2001),
                   end=c(1090, 1100, 1100, 1150, 1160, 1200, 2000))
    ir6 <- IRanges(start=c(101, 5131, 201, 301, 319, 319),
                   end=c(130, 5140, 250, 318, 318, 340))
    target <- list(chr2=ir2, chr6=ir6)
    checkIdentical(target, current)

    current <- as.list(cigarRangesAlongReferenceSpace(cigar, pos=pos, f=rname,
                                                      ops=ops,
                                                      reduce.ranges=TRUE))
    ir2b <- c(reduce(ir2[1:6]), reduce(ir2[7]))
    ir6b <- c(reduce(ir6[1:2]), reduce(ir6[3]), reduce(ir6[4:6]))
    target <- list(chr2=ir2b, chr6=ir6b)
    checkIdentical(target, current)

    current <- as.list(extractAlignmentRangesOnReference(cigar, pos=pos,
                                                         f=rname))
    checkIdentical(target, current)

    current <- as.list(cigarRangesAlongReferenceSpace(cigar, pos=pos, f=rname,
                                                      ops=setdiff(ops, "D")))
    ir2 <- ir2[-5]
    target <- list(chr2=ir2, chr6=ir6)
    checkIdentical(target, current)

    current <- as.list(cigarRangesAlongReferenceSpace(cigar, pos=pos, f=rname,
                                                      ops=setdiff(ops, "D"),
                                                      reduce.ranges=TRUE))
    ir2 <- c(reduce(ir2[1:5]), reduce(ir2[6]))
    ir6 <- c(reduce(ir6[1:2]), reduce(ir6[3]), reduce(ir6[4:6]))
    target <- list(chr2=ir2, chr6=ir6)
    checkIdentical(target, current)

    current <- as.list(extractAlignmentRangesOnReference(cigar, pos=pos,
                                                         drop.D.ranges=TRUE,
                                                         f=rname))
    checkIdentical(target, current)
}

test_cigarQNarrow <- function()
{
    cigar <- c("25M4D10M", "6S17M6I3M3S")

    ans <- cigarQNarrow(cigar)
    ans0 <- cigar
    attr(ans0, "rshift") <- c(0L, 0L)
    checkIdentical(ans, ans0)

    ans <- cigarQNarrow(cigar, start=3, end=-3)
    ans0 <- c("23M4D8M", "4S17M6I3M1S")
    attr(ans0, "rshift") <- c(2L, 0L)
    checkIdentical(ans, ans0)

    ans <- cigarQNarrow(cigar, start=7, end=-4)
    ans0 <- c("19M4D7M", "17M6I3M")
    attr(ans0, "rshift") <- c(6L, 0L)
    checkIdentical(ans, ans0)

    ans <- cigarQNarrow(cigar, start=8, end=-5)
    ans0 <- c("18M4D6M", "16M6I2M")
    attr(ans0, "rshift") <- c(7L, 1L)
    checkIdentical(ans, ans0)

    ans <- cigarQNarrow(cigar, start=25)
    ans0 <- c("1M4D10M", "5I3M3S")
    attr(ans0, "rshift") <- c(24L, 17L)
    checkIdentical(ans, ans0)

    ans <- cigarQNarrow(cigar, start=26)
    ans0 <- c("10M", "4I3M3S")
    attr(ans0, "rshift") <- c(29L, 17L)
    checkIdentical(ans, ans0)

    ans <- cigarQNarrow(cigar, start=26, end=-8)
    ans0 <- c("3M", "3I")
    attr(ans0, "rshift") <- c(29L, 17L)
    checkIdentical(ans, ans0)

    ans <- cigarQNarrow(cigar, start=26, end=-10)
    ans0 <- c("1M", "1I")
    attr(ans0, "rshift") <- c(29L, 17L)
    checkIdentical(ans, ans0)
}

test_cigarNarrow <- function()
{
    cigar <- c("25M4D10M", "6S17M6I3M3S")

    ans <- cigarNarrow(cigar)
    ans0 <- c("25M4D10M", "17M6I3M")
    attr(ans0, "rshift") <- c(0L, 0L)
    checkIdentical(ans, ans0)

    ans <- cigarNarrow(cigar, start=3, end=-3)
    ans0 <- c("23M4D8M", "15M6I1M")
    attr(ans0, "rshift") <- c(2L, 2L)
    checkIdentical(ans, ans0)

    ans <- cigarNarrow(cigar, start=7, end=-4)
    ans0 <- c("19M4D7M", "11M")
    attr(ans0, "rshift") <- c(6L, 6L)
    checkIdentical(ans, ans0)

    ans <- cigarNarrow(cigar, start=8, end=-5)
    ans0 <- c("18M4D6M", "9M")
    attr(ans0, "rshift") <- c(7L, 7L)
    checkIdentical(ans, ans0)

    ans <- cigarNarrow(cigar[1], start=26, end=-10)
    ans0 <- "1M"
    attr(ans0, "rshift") <- 29L
    checkIdentical(ans, ans0)
}

test_refLocs2queryLocs <- function() {
  cigar <- "66S42M2I20M8I18D15M43243N5M1D38M1D85M1D115M139S"
  pos <- 525842L
  ref <- 43425L + pos - 1L
  query <- 238L
  ans <- .Call("ref_locs_to_query_locs", ref, cigar, pos, FALSE,
               PACKAGE="GenomicRanges")
  checkIdentical(ans, query)
}
