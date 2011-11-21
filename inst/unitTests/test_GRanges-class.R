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

    ## elementMetadata
    checkException(elementMetadata(gr) <- DataFrame(strand = 1:length(gr)),
                   silent = TRUE)
    checkException(elementMetadata(gr) <- DataFrame(score = letters),
                   silent = TRUE)

    gr <- make_test_GRanges()
    elementMetadata(gr) <- NULL
    checkIdentical(elementMetadata(gr),
                   new("DataFrame", nrows = length(gr)))

    gr <- make_test_GRanges()
    val <- DataFrame(x = 1:length(gr), y = head(letters, length(gr)))
    elementMetadata(gr) <- val
    checkTrue(validObject(gr))
    checkIdentical(elementMetadata(gr), val)
    rownames(val) <- names(gr)
    checkIdentical(elementMetadata(gr, row.names=TRUE), val)
    elementMetadata(gr) <- val
    checkTrue(validObject(gr))
    checkIdentical(elementMetadata(gr, row.names=TRUE), val)
    rownames(val) <- NULL
    checkIdentical(elementMetadata(gr), val)

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
    checkException(resize(gr, 10, fix = "middle"), silent = TRUE)
    checkException(resize(gr, 10, fix = rep("end", 3)), silent = TRUE)
    resized <- resize(gr, 10)
    checkIdentical(rep(10L, length(gr)), width(resized))
    checkIdentical(c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 1L, 1L), start(resized))
    checkIdentical(ranges(resize(gr, 10, fix = "center")),
                   IRanges(rep(1:5, each=2), width = 10,
                           names = head(letters, 10)))
    checkIdentical(ranges(resize(gr, 10, fix = c("start", "end"))),
                   IRanges(c(1L, 1L, 3L, 1L, 5L, 1L, 7L, 1L, 1L, 10L),
                           width = 10, names = head(letters, 10)))
                   

    ## shift
    gr <- make_test_GRanges()
    shifted <- shift(gr, 10)
    checkIdentical(start(gr) + 10L, start(shifted))
    checkIdentical(width(gr), width(shifted))
    gr <- GRanges("chrA", IRanges(20, 30), seqlengths=c(chrA=100))
    checkIdentical(IRanges(8, 18), ranges(shift(gr, -12)))
    shifted <- suppressWarnings(shift(gr, 78))
    checkIdentical(IRanges(98, 100), ranges(shifted))

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

    ## restrict
    gr <-  make_test_GRanges()
    st <- structure(c(4,5), names = c("chr1", "chr2"))
    en <-  structure(c(8,9), names = c("chr2", "chr3"))
    res <- restrict(gr, start = st, end = en)
    checkIdentical(elementMetadata(gr), elementMetadata(res))
    checkIdentical(seqnames(gr), seqnames(res))
    checkIdentical(seqinfo(gr), seqinfo(res))
    target <- IRanges(start=c(4, 5, 5, 5, 5, 6, 7, 8, 9, 10),
                      end = c(10, 8, 8, 8, 10, 10, 9, 9, 9, 9),
                      names=letters[1:10])
    checkIdentical(ranges(res), target)
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

test_GRanges_precede_follow <- function() {

    g1 <- GRanges(seqnames = c("chr2", "chr2", rep("chr1", 6)),
        ranges = IRanges(c(45, 10, 10, 10, 10, 35, 35 ,35),  c(50, 15, 15, 15, 15, 40, 40, 40)),
        strand = c("+", "+", "+", "-", "*", "+", "-", "*"))

    g2 <- GRanges( seqnames = c("chr1", "chr1", "chr1", "chr2", "chr3", "chr2"),
        ranges = IRanges(c(20, 20, 20, 20, 20, 30), c( 30, 30, 30, 30, 30, 40)),
        strand = c("+", "-", "*", "+", "*", "+"))

    current <- precede(g1, g2)
    tmp <- NA_integer_
    target <-  c(tmp, 4, 1, tmp, 1, tmp ,3, 2)
    checkEquals(target, current)

    current <- follow(g1, g2)
    target  <- c(6, tmp ,tmp, 2, 2, 3, tmp, 1)
    checkEquals(target, current)
}

test_GRanges_nearest <- function() {
    g1 <- GRanges(seqnames = c("chr2", "chr2", rep("chr1", 6)),
        ranges = IRanges(c(45, 10, 10, 10, 10, 35, 35 ,35),  c(50, 15, 15, 15, 15, 40, 40, 40)),
        strand = c("+", "+", "+", "-", "*", "+", "-", "*"))

    g2 <- GRanges( seqnames = c("chr1", "chr1", "chr1", "chr2", "chr3", "chr2"),
        ranges = IRanges(c(20, 20, 20, 20, 20, 30), c( 30, 30, 30, 30, 30, 40)),
        strand = c("+", "-", "*", "+", "*", "+"))
    current <- nearest(g1, g2)
    target <- c(6, 4, 1, 2, 1, 3, 3, 3)
    checkEquals(target, current)
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
  vc1 <- as.data.frame(values(gc1))
  rownames(vc1) <- NULL
  vc.orig <- as.data.frame(rbind(values(gr1), values(gr2)))
  rownames(vc.orig) <- NULL
  checkIdentical(vc1, vc.orig)
  
  #############################################################################
  ## Combining GRanges objects with differing elementMetadata
  colnames(values(gr1))[1] <- 'illegal'
  checkException(c(gr1, gr2), silent=TRUE)
  
  ## Ignore elementMetadata
  gc2 <- c(gr1, gr2, .ignoreElementMetadata=TRUE)
  em2 <- elementMetadata(gc2)
  checkIdentical(nrow(em2), length(gc2))
  checkIdentical(ncol(em2), 0L)
}

