make_test_GRanges <- function() {
    new("GRanges",
        seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")),
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

    checkIdentical(GRanges(seqnames =
                           Rle(c("chr1", "chr2", "chr1", "chr3"),
                               c(1, 3, 2, 4)),
                           ranges = IRanges(1:10, width = 10:1,
                             names = head(letters,10)),
                           strand = Rle(factor(c("-", "+", "*", "+", "-")),
                             c(1, 2, 2, 3, 2)),
                           score = 1:10, GC = seq(1, 0, length=10)),
                   make_test_GRanges())

    seqinfo <- Seqinfo(letters, rep(1000L, length(letters)))
    checkIdentical(seqinfo(GRanges(seqinfo=seqinfo)), seqinfo)
    checkIdentical(seqinfo(GRanges(seqlengths = seqlengths(seqinfo))), seqinfo)
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

    ## mcols
    checkException(mcols(gr) <- DataFrame(strand = 1:length(gr)),
                   silent = TRUE)
    checkException(mcols(gr) <- DataFrame(score = letters),
                   silent = TRUE)

    gr <- make_test_GRanges()
    mcols(gr) <- NULL
    checkIdentical(mcols(gr),
                   new("DataFrame", nrows = length(gr)))

    gr <- make_test_GRanges()
    val <- DataFrame(x = 1:length(gr), y = head(letters, length(gr)))
    mcols(gr) <- val
    checkTrue(validObject(gr))
    checkIdentical(mcols(gr), val)
    rownames(val) <- names(gr)
    checkIdentical(mcols(gr, row.names=TRUE), val)
    mcols(gr) <- val
    checkTrue(validObject(gr))
    checkIdentical(mcols(gr, row.names=TRUE), val)
    rownames(val) <- NULL
    checkIdentical(mcols(gr), val)

    ## names
    checkException(names(gr) <- letters, silent = TRUE)

    gr <- make_test_GRanges()
    names(gr) <- NULL
    checkIdentical(names(gr), NULL)

    gr <- make_test_GRanges()
    names(gr) <- head(letters, length(gr))
    checkIdentical(names(gr), head(letters, length(gr)))

    ## seqlevels
    gr <- make_test_GRanges()
    val <- seqlevels(gr)
    val <- gsub("chr", "Chr", val)
    seqlevels(gr) <- val
    checkIdentical(seqlevels(gr), val)

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
}

test_GRanges_Vector <- function() {
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
    checkIdentical(gr, gr[Rle(TRUE)])

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
    checkIdentical(as.data.frame(c(gr, gr2), row.names=NULL),
                   rbind(as.data.frame(gr, row.names=NULL), as.data.frame(gr2, row.names=NULL)))
    checkIdentical(as.data.frame(c(gr2, gr), row.names=NULL),
                   rbind(as.data.frame(gr2, row.names=NULL), as.data.frame(gr, row.names=NULL)))

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
    checkIdentical(split(gr, rep(c("a", "b"), each=5)),
                   GRangesList(a = head(gr, 5), b = tail(gr, 5)))

    ## window
    gr <- make_test_GRanges()
    checkIdentical(gr[1:3], window(gr, 1, 3))
}

test_GRanges_combine <- function() {
  gr1 <- make_test_GRanges()
  gr2 <- make_test_GRanges()
  
  #############################################################################
  ## An unremarkable combination
  gc1 <- c(gr1, gr2)
  checkEquals(start(gc1), c(start(gr1), start(gr2)))
  checkEquals(end(gc1), c(end(gr1), end(gr2)))
  
  ## Check the combined data frames -- the rownaming is different when
  ## combining using these two strategies, so ignore them for now.
  vc1 <- as.data.frame(mcols(gc1))
  rownames(vc1) <- NULL
  vc.orig <- as.data.frame(rbind(mcols(gr1), mcols(gr2)))
  rownames(vc.orig) <- NULL
  checkIdentical(vc1, vc.orig)
  
  #############################################################################
  ## Combining GRanges objects with differing metadata columns
  colnames(mcols(gr1))[1] <- 'illegal'
  checkException(c(gr1, gr2), silent=TRUE)
  
  ## Ignore mcols
  gc2 <- c(gr1, gr2, ignore.mcols=TRUE)
  em2 <- mcols(gc2)
  checkIdentical(nrow(em2), length(gc2))
  checkIdentical(ncol(em2), 0L)
}

