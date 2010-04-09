make_test_GRanges <- function() {
    new("GRanges",
        seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        seqlengths =
        c("chr1" = NA_integer_, "chr2" = NA_integer_, "chr3" = NA_integer_),
        elementMetadata = DataFrame(score = 1:10, GC = seq(1, 0, length=10)))
}

test_GRanges_construction <- function() {
    checkException(GRanges(letters), silent = TRUE)
    checkException(GRanges(ranges = IRanges(1:10, 1:10)), silent = TRUE)
    checkException(GRanges(letters, IRanges(1:10, 1:10)), silent = TRUE)
    checkException(GRanges(letters, IRanges(1:26, 1:26), strand = letters),
                   silent = TRUE)
    checkException(GRanges(letters, IRanges(1:26, 1:26), score = 1:10),
                   silent = TRUE)
    checkException(GRanges(letters, IRanges(1:26, 1:26), start = 1:26),
                   silent = TRUE)
    checkException(GRanges(letters, IRanges(1:26, 1:26), end = 1:26),
                   silent = TRUE)
    checkException(GRanges(letters, IRanges(1:26, 1:26), width = 1:26),
                   silent = TRUE)
    checkException(GRanges(letters, IRanges(1:26, 1:26), element = letters),
                   silent = TRUE)
    checkException(GRanges(c(letters, NA), IRanges(1:27, 1:27)),
                   silent = TRUE)

    checkTrue(validObject(new("GRanges")))
    checkTrue(validObject(GRanges()))
    checkTrue(validObject(GRanges(letters, IRanges(1:26, 1:26))))
    checkTrue(validObject(GRanges(letters, IRanges(1:26, 1:26), score = 1:26)))
    checkTrue(validObject(GRanges(factor(letters), IRanges(1:26, 1:26))))
    checkTrue(validObject(GRanges(1:10, IRanges(1:10, 1:10))))

    checkIdentical(GRanges(seqnames =
                           Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
                           ranges =
                           IRanges(1:10, width = 10:1, names = head(letters,10)),
                           strand =
                           Rle(strand(c("-", "+", "*", "+", "-")),
                               c(1, 2, 2, 3, 2)),
                           score = 1:10, GC = seq(1, 0, length=10)),
                   make_test_GRanges())
}

test_GRanges_coercion <- function() {
    ## no strand or score
    gr <-
      GRanges(seqnames = c(1,1,2),
              ranges = IRanges(1:3,4:6, names = head(letters,3)))
    df <-
      data.frame(seqnames = factor(c(1,1,2)),
                 start = 1:3, end = 4:6, width = c(4L, 4L, 4L),
                 strand = strand(rep("*", 3)),
                 row.names = head(letters,3),
                 stringsAsFactors = FALSE)
    checkIdentical(as.data.frame(gr), df)

    ## score, no strand
    gr <-
      GRanges(seqnames = c(1,1,2),
              ranges = IRanges(1:3,4:6, names = head(letters,3)),
              score = c(10L,2L,NA))
    df <-
      data.frame(seqnames = factor(c(1,1,2)),
                 start = 1:3, end = 4:6, width = c(4L, 4L, 4L),
                 strand = strand(rep("*", 3)),
                 score = c(10L,2L,NA),
                 row.names = head(letters,3),
                 stringsAsFactors = FALSE)
    checkIdentical(as.data.frame(gr), df)

    ## strand, no score
    gr <-
      GRanges(seqnames = c(1,1,2),
              ranges = IRanges(1:3,4:6, names = head(letters,3)),
              strand = strand(c("+", "-", "*")))
    df <-
      data.frame(seqnames = factor(c(1,1,2)),
                 start = 1:3, end = 4:6, width = c(4L, 4L, 4L),
                 strand = strand(c("+", "-", "*")),
                 row.names = head(letters,3),
                 stringsAsFactors = FALSE)
    checkIdentical(as.data.frame(gr), df)

    ## strand & score
    gr <-
      GRanges(seqnames = c(1,1,2),
              ranges = IRanges(1:3,4:6, names = head(letters,3)),
              strand = strand(c("+", "-", "*")),
              score = c(10L,2L,NA))
    df <-
      data.frame(seqnames = factor(c(1,1,2)),
                 start = 1:3, end = 4:6, width = c(4L, 4L, 4L),
                 strand = strand(c("+", "-", "*")),
                 score = c(10L,2L,NA),
                 row.names = head(letters,3),
                 stringsAsFactors = FALSE)
    checkIdentical(as.data.frame(gr), df)
}

test_GRanges_accessors <- function() {
    ## seqnames
    checkException(seqnames(GRanges()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GRanges()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GRanges()) <- letters,
                   silent = TRUE)

    gr <- make_test_GRanges()
    val <- seqnames(gr)
    runValue(val) <- factor(paste(runValue(val), ".new", sep=""))
    seqnames(gr) <- val
    checkIdentical(seqnames(gr), val)

    gr <- make_test_GRanges()
    val <- gsub("chr", "Chr", seqnames(gr))
    seqnames(gr) <- val
    checkIdentical(seqnames(gr), val)

    ## ranges
    checkException(ranges(GRanges()) <- NULL, silent = TRUE)
    checkException(ranges(make_test_GRanges()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GRanges()) <- IRanges(1:26, 1:26),
                   silent = TRUE)

    gr <- make_test_GRanges()
    val <- IRanges(1:length(gr), width = 10)
    ranges(gr) <- val
    checkIdentical(ranges(gr), val)

    ## strand
    checkException(strand(GRanges()) <- NULL, silent = TRUE)
    checkException(strand(make_test_GRanges()) <- NULL, silent = TRUE)
    checkException(strand(make_test_GRanges()) <- letters, silent = TRUE)

    gr <- make_test_GRanges()
    val <- Rle(strand("+"), length(gr))
    strand(gr) <- val
    checkIdentical(strand(gr), val)

    gr <- make_test_GRanges()
    val <- rep(strand("+"), length(gr))
    strand(gr) <- val
    checkIdentical(strand(gr), Rle(val))

    ## seqlengths
    checkException(seqlengths(GRanges()) <- NULL, silent = TRUE)
    checkException(seqlengths(make_test_GRanges()) <- NULL, silent = TRUE)
    checkException(seqlengths(make_test_GRanges()) <- 1:10,
                   silent = TRUE)

    gr <- make_test_GRanges()
    val <- seqlengths(gr)
    val[] <- c(10L, 20L, 30L)
    seqlengths(gr) <- val
    checkIdentical(seqlengths(gr), val)

    ## elementMetadata
    checkException(elementMetadata(gr) <- DataFrame(strand = 1:length(gr)),
                   silent = TRUE)
    checkException(elementMetadata(gr) <- DataFrame(score = letters),
                   silent = TRUE)

    gr <- make_test_GRanges()
    elementMetadata(gr) <- NULL
    checkIdentical(elementMetadata(gr),
                   new("DataFrame", nrows = length(gr), rownames = names(gr)))

    gr <- make_test_GRanges()
    val <- DataFrame(x = 1:length(gr), y = head(letters, length(gr)))
    rownames(val) <- names(gr)
    elementMetadata(gr) <- val
    checkTrue(validObject(gr))
    checkIdentical(elementMetadata(gr), val)

    ## names
    checkException(names(gr) <- letters, silent = TRUE)

    gr <- make_test_GRanges()
    names(gr) <- NULL
    checkIdentical(names(gr), NULL)

    gr <- make_test_GRanges()
    names(gr) <- head(letters, length(gr))
    checkIdentical(names(gr), head(letters, length(gr)))
}

test_GRanges_Ranges <- function() {
    ## start
    checkException(start(GRanges()) <- NULL, silent = TRUE)
    checkException(start(make_test_GRanges()) <- letters, silent = TRUE)
    checkException(start(make_test_GRanges()) <- 1:26, silent = TRUE)

    gr <- make_test_GRanges()
    start(gr) <- as.numeric(seq_len(length(gr)))
    checkIdentical(start(gr), seq_len(length(gr)))

    ## end
    checkException(end(GRanges()) <- NULL, silent = TRUE)
    checkException(end(make_test_GRanges()) <- letters, silent = TRUE)
    checkException(end(make_test_GRanges()) <- 1:26, silent = TRUE)

    gr <- make_test_GRanges()
    end(gr) <- as.numeric(10L + seq_len(length(gr)))
    checkIdentical(end(gr), 10L + seq_len(length(gr)))

    ## width
    checkException(width(GRanges()) <- NULL, silent = TRUE)
    checkException(width(make_test_GRanges()) <- letters, silent = TRUE)
    checkException(width(make_test_GRanges()) <- 1:26, silent = TRUE)

    gr <- make_test_GRanges()
    width(gr) <- as.numeric(10L + seq_len(length(gr)))
    checkIdentical(width(gr), 10L + seq_len(length(gr)))

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
    checkException(resize(gr, 10, fix = "end"), silent = TRUE)
    resized <- resize(gr, 10)
    checkIdentical(rep(10L, length(gr)), width(resized))
    checkIdentical(c(1L, 2L, 3L, 2L, 3L, 6L, 7L, 8L, 1L, 1L), start(resized))

    ## shift
    gr <- make_test_GRanges()
    shifted <- shift(gr, 10)
    checkIdentical(start(gr) + 10L, start(shifted))
    checkIdentical(width(gr), width(shifted))

    ## disjoin
    gr <- unname(make_test_GRanges())[ , character(0)]
    checkIdentical(disjoin(gr),
                   GRanges(seqnames = Rle(c("chr1", "chr2", "chr3"), c(3, 3, 4)),
                           ranges = IRanges(start=c(6, 1, 5, 2, 3, 4, 7, 8, 9, 10),
                                            end=c(10, 10, 10, 2, 10, 10, 7, 10, 9, 10)),
                           strand =
                           strand(c("+", "-", "*", "+", "+", "*", "+", "+", "-", "-"))))

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

    ## coverage
    gr <- make_test_GRanges()
    checkIdentical(coverage(gr),
                   RleList("chr1" = Rle(1:3, c(4, 1, 5)),
                           "chr2" = Rle(0:3, c(1, 1, 1, 7)),
                           "chr3" = Rle(0:4, c(6, 1, 1, 1, 1))))
    checkIdentical(coverage(gr, width = list(10, 20, 30)),
                   RleList("chr1" = Rle(1:3, c(4, 1, 5)),
                           "chr2" = Rle(c(0:3, 0L), c(1, 1, 1, 7, 10)),
                           "chr3" = Rle(c(0:4, 0L), c(6, 1, 1, 1, 1, 20))))
    checkIdentical(coverage(gr, weight = list(1L, 10L, 100L)),
                   RleList("chr1" = Rle(1:3, c(4, 1, 5)),
                           "chr2" = Rle(10L * 0:3, c(1, 1, 1, 7)),
                           "chr3" = Rle(100L * 0:4, c(6, 1, 1, 1, 1))))
    checkIdentical(coverage(gr, shift = list(0, 1, 2)),
                   RleList("chr1" = Rle(1:3, c(4, 1, 5)),
                           "chr2" = Rle(0:3, c(2, 1, 1, 7)),
                           "chr3" = Rle(0:4, c(8, 1, 1, 1, 1))))
}

test_GRanges_Sequence <- function() {
    ## [
    gr <- make_test_GRanges()
    checkException(gr[1000], silent = TRUE)
    checkException(gr["bad"], silent = TRUE)
    checkIdentical(gr, gr[])
    checkIdentical(as.data.frame(gr)[c(1,3,5),], as.data.frame(gr[c(1,3,5)]))
    checkIdentical(as.data.frame(gr)[c(1,3,5),-7],
                   as.data.frame(gr[c(1,3,5),"score"]))
    checkIdentical(as.data.frame(gr)[c(1,3,5),-7],
                   as.data.frame(gr[c(1,3,5),1]))

    ## [<-
    gr <- make_test_GRanges()
    gr[] <- rev(gr)
    revgr <- rev(make_test_GRanges())
    names(revgr) <- rev(names(revgr))
    checkIdentical(gr, revgr)

    ## c
    gr <- make_test_GRanges()
    gr2 <- gr
    names(gr2) <- NULL
    checkException(c(GRanges(), RangedData()), silent = TRUE)
    checkException(c(gr, gr[,-1]), silent = TRUE)
    checkIdentical(as.data.frame(c(gr, gr)),
                   as.data.frame(gr)[rep(seq_len(length(gr)), 2),])
    checkIdentical(as.data.frame(c(gr, gr2)),
                   rbind(as.data.frame(gr), as.data.frame(gr2)))
    checkIdentical(as.data.frame(c(gr2, gr)),
                   rbind(as.data.frame(gr2), as.data.frame(gr)))

    ## length
    checkIdentical(length(gr), length(gr@seqnames))

    ## seqselect
    gr <- make_test_GRanges()
    checkIdentical(gr[1:3], seqselect(gr, 1, 3))
    checkIdentical(gr[c(1:3, 1:3)], seqselect(gr, c(1,1), c(3,3)))

    ## seqselect<-
    gr1 <- make_test_GRanges()
    gr1[1:3] <- make_test_GRanges()[4:6]
    gr2 <- make_test_GRanges()
    seqselect(gr2, 1, 3) <- make_test_GRanges()[4:6]
    checkIdentical(gr1, gr2)

    ## split
    gr <- make_test_GRanges()
    checkException(split(gr, NULL), silent = TRUE)
    checkIdentical(split(unname(gr)),
                   GRangesList(lapply(structure(seq_len(length(gr)),
                                names = as.character(seq_len(length(gr)))),
                                      function(i) unname(gr)[i])))
    checkIdentical(split(gr),
                   GRangesList(lapply(structure(seq_len(length(gr)),
                                                names = names(gr)),
                                      function(i) gr[i])))
    checkIdentical(split(gr, rep(c("a", "b"), each=5)),
                   GRangesList(a = head(gr, 5), b = tail(gr, 5)))

    ## window
    gr <- make_test_GRanges()
    checkIdentical(gr[1:3], window(gr, 1, 3))
}
