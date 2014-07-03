###

.TARGET_names <- letters[1:10]
.TARGET_seqlevels <- c("chr1", "chr2", "chr3", "chrX")
.TARGET_seqnames <- Rle(factor(c("chr1", "chr2", "chr1", "chr3"),
                               levels=.TARGET_seqlevels),
                        c(1, 3, 2, 4))
.TARGET_start <- 1:10
.TARGET_end <- rep(10L, 10)
.TARGET_ranges <- IRanges(.TARGET_start, .TARGET_end, names=.TARGET_names)
.TARGET_strand <- Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2))
.TARGET_mcols <- DataFrame(score=1:10, GC=seq(1, 0, length=10))
.TARGET_seqlengths <- c(120L, NA, 80L, 50L)
names(.TARGET_seqlengths) <- .TARGET_seqlevels
.TARGET_seqinfo <- Seqinfo(seqnames=.TARGET_seqlevels,
                           seqlengths=.TARGET_seqlengths)

.make_TARGET_GRanges <- function()
{
    new("GRanges", seqnames=.TARGET_seqnames,
                   ranges=.TARGET_ranges,
                   strand= .TARGET_strand,
                   elementMetadata=.TARGET_mcols,
                   seqinfo=.TARGET_seqinfo)
}

test_GRanges_construction <- function()
{
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

    current_seqnames <- S4Vectors:::decodeRle(.TARGET_seqnames)
    current_strand <- Rle(c("-", "+", "*", "+", "-"),
                          c(1, 2, 2, 3, 2))
    current <- GRanges(seqnames=current_seqnames,
                       ranges=.TARGET_ranges,
                       strand=current_strand,
                       seqlengths=.TARGET_seqlengths,
                       score=1:10, GC=seq(1, 0, length=10))
    checkIdentical(.make_TARGET_GRanges(), current)

    ## Call with unnamed 'seqnames', 'ranges', and 'strand' args.
    current <- GRanges(current_seqnames,
                       .TARGET_ranges,
                       current_strand,
                       seqlengths=.TARGET_seqlengths,
                       score=1:10, GC=seq(1, 0, length=10))
    checkIdentical(.make_TARGET_GRanges(), current)

    ## Call with unnamed metadata cols.
    score <- 1:10
    GC <- seq(1, 0, length=10)
    current <- GRanges(.TARGET_seqnames,
                       .TARGET_ranges,
                       .TARGET_strand,
                       seqlengths=.TARGET_seqlengths,
                       score, GC)
    checkIdentical(.make_TARGET_GRanges(), current)

    ## Call with 'c' metadata col. 
    current <- GRanges(seqnames=.TARGET_seqnames,
                       ranges=.TARGET_ranges,
                       strand=.TARGET_strand,
                       c=LETTERS[1:10])
    checkIdentical(DataFrame(c=LETTERS[1:10]), mcols(current))

    seqinfo <- Seqinfo(letters, rep(1000L, length(letters)))
    checkIdentical(seqinfo(GRanges(seqinfo=seqinfo)), seqinfo)
    checkIdentical(seqinfo(GRanges(seqlengths = seqlengths(seqinfo))), seqinfo)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors.
###

test_GRanges_length <- function()
{
    checkIdentical(length(.TARGET_ranges), length(.make_TARGET_GRanges()))
}

test_GRanges_names <- function()
{
    ## names() getter
    gr0 <- GRanges()
    gr1 <- .make_TARGET_GRanges()
    checkIdentical(NULL, names(gr0))
    checkIdentical(.TARGET_names, names(gr1))

    ## names() setter
    names(gr0) <- names(gr0)  # no-op
    checkIdentical(GRanges(), gr0)
    names(gr1) <- names(gr1)  # no-op
    checkIdentical(.make_TARGET_GRanges(), gr1)
    checkException(names(gr1) <- letters, silent = TRUE)
    names(gr1) <- NULL
    checkIdentical(NULL, names(gr1))
}

test_GRanges_seqnames <- function()
{
    ## seqnames() getter
    gr0 <- GRanges()
    gr1 <- .make_TARGET_GRanges()
    checkIdentical(Rle(factor()), seqnames(gr0))
    checkIdentical(.TARGET_seqnames, seqnames(gr1))

    ## seqnames() setter
    seqnames(gr0) <- seqnames(gr0)  # no-op
    checkIdentical(GRanges(), gr0)
    seqnames(gr1) <- seqnames(gr1)  # no-op
    checkIdentical(.make_TARGET_GRanges(), gr1)
    checkException(seqnames(gr0) <- NULL, silent=TRUE)
    checkException(seqnames(gr1) <- NULL, silent=TRUE)
    checkException(seqnames(gr1) <- letters, silent=TRUE)
}

test_GRanges_ranges <- function()
{
    ## ranges() getter
    gr0 <- GRanges()
    gr1 <- .make_TARGET_GRanges()
    checkIdentical(IRanges(), ranges(gr0))
    checkIdentical(.TARGET_ranges, ranges(gr1))

    ## ranges() setter
    ranges(gr0) <- ranges(gr0)  # no-op
    checkIdentical(GRanges(), gr0)
    ranges(gr1) <- ranges(gr1)  # no-op
    checkIdentical(.make_TARGET_GRanges(), gr1)
    checkException(ranges(gr0) <- NULL, silent=TRUE)
    checkException(ranges(gr1) <- NULL, silent=TRUE)
    checkException(ranges(gr1) <- IRanges(1:26, 1:26), silent=TRUE)

    val <- IRanges(1:length(gr1), width=10)
    ranges(gr1) <- val
    checkIdentical(ranges(gr1), val)
}

test_GRanges_start <- function()
{
    ## start() getter
    gr0 <- GRanges()
    gr1 <- .make_TARGET_GRanges()
    checkIdentical(integer(), start(gr0))
    checkIdentical(.TARGET_start, start(gr1))

    ## start() setter
    start(gr0) <- start(gr0)  # no-op
    checkIdentical(GRanges(), gr0)
    start(gr1) <- start(gr1)  # no-op
    checkIdentical(.make_TARGET_GRanges(), gr1)
    #checkException(start(gr0) <- NULL, silent=TRUE)    this actually works!
    checkException(suppressWarnings(start(gr1) <- letters), silent=TRUE)
    #checkException(start(gr1) <- 1:26, silent=TRUE)    this actually works!

    start(gr1) <- as.numeric(seq_len(length(gr1)))
    checkIdentical(seq_len(length(gr1)), start(gr1))
}

test_GRanges_end <- function()
{
    ## end() getter
    gr0 <- GRanges()
    gr1 <- .make_TARGET_GRanges()
    checkIdentical(integer(), end(gr0))
    checkIdentical(.TARGET_end, end(gr1))

    ## end() setter
    end(gr0) <- end(gr0)  # no-op
    checkIdentical(GRanges(), gr0)
    end(gr1) <- end(gr1)  # no-op
    checkIdentical(.make_TARGET_GRanges(), gr1)
    #checkException(end(gr0) <- NULL, silent=TRUE)    this actually works!
    checkException(suppressWarnings(end(gr1) <- letters), silent=TRUE)
    #checkException(end(gr1) <- 1:26, silent=TRUE)    this actually works!

    end(gr1) <- as.numeric(10L + seq_len(length(gr1)))
    checkIdentical(10L + seq_len(length(gr1)), end(gr1))
}

test_GRanges_width <- function()
{
    ## width() getter
    gr0 <- GRanges()
    gr1 <- .make_TARGET_GRanges()
    checkIdentical(integer(), width(gr0))
    checkIdentical(.TARGET_end - .TARGET_start + 1L, width(gr1))

    ## width() setter
    width(gr0) <- width(gr0)  # no-op
    checkIdentical(GRanges(), gr0)
    width(gr1) <- width(gr1)  # no-op
    checkIdentical(.make_TARGET_GRanges(), gr1)
    #checkException(width(gr0) <- NULL, silent=TRUE)    this actually works!
    checkException(suppressWarnings(width(gr1) <- letters), silent=TRUE)
    #checkException(width(gr1) <- 1:26, silent=TRUE)    this actually works!

    width(gr1) <- as.numeric(10L + seq_len(length(gr1)))
    checkIdentical(10L + seq_len(length(gr1)), width(gr1))
}

test_GRanges_strand <- function()
{
    ## strand() getter
    gr0 <- GRanges()
    gr1 <- .make_TARGET_GRanges()
    checkIdentical(Rle(strand()), strand(gr0))
    checkIdentical(.TARGET_strand, strand(gr1))

    ## strand() setter
    strand(gr0) <- strand(gr0)  # no-op
    checkIdentical(GRanges(), gr0)
    strand(gr1) <- strand(gr1)  # no-op
    checkIdentical(.make_TARGET_GRanges(), gr1)
    checkException(strand(gr0) <- NULL, silent=TRUE)
    checkException(strand(gr1) <- NULL, silent=TRUE)
    checkException(strand(gr1) <- letters, silent=TRUE)

    val <- Rle(strand("+"), length(gr1))
    strand(gr1) <- val
    checkIdentical(val, strand(gr1))

    strand(gr1) <- "*"
    checkIdentical(Rle(strand("*"), length(gr1)), strand(gr1))
}

test_GRanges_mcols <- function()
{
    ## mcols() getter
    gr0 <- GRanges()
    gr1 <- .make_TARGET_GRanges()
    checkIdentical(DataFrame(), mcols(gr0))
    checkIdentical(.TARGET_mcols, mcols(gr1))

    ## mcols() setter
    mcols(gr0) <- mcols(gr0)  # no-op
    checkIdentical(GRanges(), gr0)
    mcols(gr1) <- mcols(gr1)  # no-op
    checkIdentical(.make_TARGET_GRanges(), gr1)
    checkException(mcols(gr1) <- DataFrame(strand=1:length(gr)),
                   silent=TRUE)
    checkException(mcols(gr1) <- DataFrame(score=letters),
                   silent=TRUE)

    mcols(gr1) <- NULL
    checkIdentical(new("DataFrame", nrows=length(gr1)), mcols(gr1))

    val <- DataFrame(x=1:length(gr1), y = head(letters, length(gr1)))
    mcols(gr1) <- val
    checkTrue(validObject(gr1))
    checkIdentical(val, mcols(gr1))

    rownames(val) <- names(gr1)
    checkIdentical(val, mcols(gr1, use.names=TRUE))

    mcols(gr1) <- val
    checkTrue(validObject(gr1))
    checkIdentical(val, mcols(gr1, use.names=TRUE))
    rownames(val) <- NULL
    checkIdentical(val, mcols(gr1))
}

test_GRanges_seqlevels <- function()
{
    ## seqlevels() getter
    gr0 <- GRanges()
    gr1 <- .make_TARGET_GRanges()
    checkIdentical(character(), seqlevels(gr0))
    checkIdentical(.TARGET_seqlevels, seqlevels(gr1))

    ## seqlevels() setter
    seqlevels(gr0) <- seqlevels(gr0)  # no-op
    checkIdentical(GRanges(), gr0)
    seqlevels(gr1) <- seqlevels(gr1)  # no-op
    checkIdentical(.make_TARGET_GRanges(), gr1)

    val <- seqlevels(gr1)
    val <- sub("^chr", "Chr", val)
    seqlevels(gr1) <- val
    checkIdentical(val, seqlevels(gr1))
}

test_GRanges_seqlengths <- function()
{
    ## seqlengths() getter
    gr0 <- GRanges()
    gr1 <- .make_TARGET_GRanges()
    checkIdentical(setNames(integer(), character()), seqlengths(gr0))
    checkIdentical(.TARGET_seqlengths, seqlengths(gr1))

    ## seqlengths() setter
    seqlengths(gr0) <- seqlengths(gr0)  # no-op
    checkIdentical(GRanges(), gr0)
    seqlengths(gr1) <- seqlengths(gr1)  # no-op
    checkIdentical(.make_TARGET_GRanges(), gr1)
    checkException(seqlengths(gr0) <- NULL, silent=TRUE)
    checkException(seqlengths(gr1) <- NULL, silent=TRUE)
    checkException(seqlengths(gr1) <- 1:10, silent=TRUE)

    val <- seqlengths(gr1)
    val[] <- c(10L, 20L, 30L, 10L)
    seqlengths(gr1) <- val
    checkIdentical(val, seqlengths(gr1))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

test_GRanges_coercion <- function()
{
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

test_GRanges_subsetting <- function()
{
    ## [
    gr <- .make_TARGET_GRanges()
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
    gr <- .make_TARGET_GRanges()
    gr[] <- rev(gr)
    revgr <- rev(.make_TARGET_GRanges())
    names(revgr) <- rev(names(revgr))
    checkIdentical(gr, revgr)

    ## window
    gr <- .make_TARGET_GRanges()
    checkIdentical(gr[1:3], window(gr, 1, 3))

    ## [ by IRanges
    gr <- .make_TARGET_GRanges()
    checkIdentical(gr[1:3], gr[IRanges(1, 3)])
    checkIdentical(gr[c(1:3, 1:3)], gr[IRanges(c(1,1), c(3,3))])

    ## [<- by IRanges
    gr1 <- .make_TARGET_GRanges()
    gr1[1:3] <- .make_TARGET_GRanges()[4:6]
    gr2 <- .make_TARGET_GRanges()
    gr2[IRanges(1, 3)] <- .make_TARGET_GRanges()[4:6]
    checkIdentical(gr1, gr2)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine and split.
###

test_GRanges_combine <- function()
{
    gr1 <- .make_TARGET_GRanges()
    gr2 <- .make_TARGET_GRanges()
  
    #########################################################################
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

    #########################################################################
    ## Combining GRanges objects with differing metadata columns
    colnames(mcols(gr1))[1] <- 'illegal'
    checkException(c(gr1, gr2), silent=TRUE)
  
    ## Ignore mcols
    gc2 <- c(gr1, gr2, ignore.mcols=TRUE)
    em2 <- mcols(gc2)
    checkIdentical(nrow(em2), length(gc2))
    checkIdentical(ncol(em2), 0L)

    #########################################################################
    ## More testing
    gr <- .make_TARGET_GRanges()
    gr2 <- gr
    names(gr2) <- NULL
    checkException(c(GRanges(), RangedData()), silent = TRUE)
    checkException(c(gr, gr[,-1]), silent = TRUE)
    checkIdentical(as.data.frame(c(gr, gr2), row.names=NULL),
                   rbind(as.data.frame(gr, row.names=NULL),
                         as.data.frame(gr2, row.names=NULL)))
    checkIdentical(as.data.frame(c(gr2, gr), row.names=NULL),
                   rbind(as.data.frame(gr2, row.names=NULL),
                         as.data.frame(gr, row.names=NULL)))
}

test_GRanges_split <- function()
{
    gr <- .make_TARGET_GRanges()
    checkException(split(gr, NULL), silent = TRUE)
    checkIdentical(split(gr, rep(c("a", "b"), each=5)),
                   GRangesList(a = head(gr, 5), b = tail(gr, 5)))
}

