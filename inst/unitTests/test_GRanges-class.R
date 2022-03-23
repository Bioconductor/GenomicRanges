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
    checkIdentical(DataFrame(c=LETTERS[1:10]), mcols(current, use.names=FALSE))

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
    checkIdentical(.TARGET_mcols, mcols(gr1, use.names=FALSE))
    checkIdentical(`rownames<-`(.TARGET_mcols, .TARGET_names), mcols(gr1))

    ## mcols() setter
    mcols(gr0) <- mcols(gr0)  # no-op
    checkIdentical(GRanges(), gr0)
    mcols(gr1) <- mcols(gr1)  # no-op
    checkIdentical(.make_TARGET_GRanges(), gr1)
    checkException(mcols(gr1) <- DataFrame(score=letters),
                   silent=TRUE)

    mcols(gr1) <- NULL
    target <- make_zero_col_DFrame(length(gr1))
    checkIdentical(target, mcols(gr1, use.names=FALSE))
    checkIdentical(`rownames<-`(target, .TARGET_names), mcols(gr1))

    val <- DataFrame(x=1:length(gr1), y = head(letters, length(gr1)))
    mcols(gr1) <- val
    checkTrue(validObject(gr1))
    checkIdentical(val, mcols(gr1, use.names=FALSE))

    rownames(val) <- names(gr1)
    checkIdentical(val, mcols(gr1, use.names=TRUE))

    mcols(gr1) <- val
    checkTrue(validObject(gr1))
    checkIdentical(val, mcols(gr1, use.names=TRUE))
    rownames(val) <- NULL
    checkIdentical(val, mcols(gr1, use.names=FALSE))
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
  ## -- From GRanges to character/factor/data.frame -- ##

    ## no strand, no score
    gr <-
      GRanges(seqnames = factor(c("chr2", "chr11", "chr1"),
                                levels=c("chr1", "chr2", "chr11")),
              ranges = IRanges(1:-1, c(4:5, -2), names=head(letters, 3)))
    target1 <- c(a="chr2:1-4", b="chr11:0-5", c="chr1:-1--2")
    checkIdentical(target1, as.character(gr))
    target2 <- factor(target1, levels=target1[c(3, 1, 2)])
    checkIdentical(target2, as.factor(gr))
    target3 <-
      data.frame(seqnames = factor(c("chr2", "chr11", "chr1"),
                                   levels=c("chr1", "chr2", "chr11")),
                 start = 1:-1, end = c(4:5, -2L), width = c(4L, 6L, 0L),
                 strand = strand(rep("*", 3)),
                 row.names = head(letters, 3),
                 stringsAsFactors = FALSE)
    checkIdentical(target3, as.data.frame(gr))

    ## strand, no score
    strand(gr) <- strand(c("+", "-", "*"))
    target1s <- c(a="chr2:1-4:+", b="chr11:0-5:-", c="chr1:-1--2:*")
    checkIdentical(target1s, as.character(gr))
    checkIdentical(target1, as.character(gr, ignore.strand=TRUE))
    target2s <- factor(target1s, levels=target1s[c(3, 1, 2)])
    checkIdentical(target2s, as.factor(gr))
    target3$strand <- strand(c("+", "-", "*"))
    checkIdentical(target3, as.data.frame(gr))

    ## strand, score
    mcols(gr)$score <- c(10L, 2L, NA)
    checkIdentical(target1s, as.character(gr))
    checkIdentical(target2s, as.factor(gr))
    target3$score <- c(10L, 2L, NA)
    checkIdentical(target3, as.data.frame(gr))

    ## no strand, score
    gr <- unstrand(gr)
    checkIdentical(target1, as.character(gr))
    checkIdentical(target2, as.factor(gr))
    target3$strand <- strand("*")
    checkIdentical(target3, as.data.frame(gr))

  ## -- From GRanges to character/factor (continued) -- ##

    checkIdentical(as.character(gr), as(gr, "character"))
    checkIdentical(as.factor(gr), as(gr, "factor"))

    set.seed(555)
    gr2 <- sample(gr, 100, replace=TRUE)
    current <- as.factor(gr2)
    checkTrue(all(as.factor(as.character(gr2)) == current))
    checkIdentical(unname(as.character(sort(unique(gr2)))), levels(current))

    strand(gr2) <- c("*", "-", "+", "-")
    current <- as.factor(gr2)
    checkTrue(all(as.factor(as.character(gr2)) == current))
    checkIdentical(unname(as.character(sort(unique(gr2)))), levels(current))

  ## -- table() -- ##

    current <- table(gr2)  # same as table(as.factor(gr2)) but much faster
    target <- table(as.factor(gr2))
    dimnames(current) <- unname(dimnames(current))
    dimnames(target) <- unname(dimnames(target))
    checkIdentical(target, current)

  ## -- From character/factor to GRanges -- ##

    x <- c(a="chrX:21-3.5e+03",
           b="chr1: +21 \t-+30.5:-",
           c="1:-15--3:*",
           d="chr..Y:-21..-3:",
           e="chr1-X: \t 21 ..  \t+30  \t\t:+")
    current <- as(x, "GRanges")
    target <- GRanges(c("chrX", "chr1", "1", "chr..Y", "chr1-X"),
                      IRanges(c(  21, 21, -15, -21, 21),
                              c(3500, 30,  -3,  -3, 30),
                              names=letters[1:5]),
                      c("*", "-", "*", "*", "+"))
    checkIdentical(target, current)
    checkIdentical(target, as(x, "GenomicRanges"))
    checkIdentical(unname(target), as(unname(x), "GRanges"))
    f <- as.factor(x)
    checkIdentical(target, as(f, "GRanges"))
    checkIdentical(target, as(f, "GenomicRanges"))
    checkIdentical(unname(target), as(unname(f), "GRanges"))

  ## -- Going back and forth between character/factor and GRanges -- ##

    ## this looses the metadata(), mcols(), and seqinfo()
    current <- as(as.character(gr2), "GRanges")
    target <- gr2
    mcols(target) <- NULL
    seqlevels(target) <- seqlevels(current)
    checkIdentical(target, current)
    checkIdentical(target, as(as.factor(gr2), "GRanges"))

    current <- as(as.character(gr2, ignore.strand=TRUE), "GRanges")
    checkIdentical(unstrand(target), current)

    ##
    x <- as.character(gr2)
    checkIdentical(x, as.character(as(x, "GRanges")))
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
### Concatenation
###

test_GRanges_concatenate <- function()
{
    gr1 <- gr2 <- .make_TARGET_GRanges()

    #########################################################################
    ## An unremarkable concatenation
    gr12 <- c(gr1, gr2)
    checkIdentical(seqnames(gr12), c(seqnames(gr1), seqnames(gr2)))
    checkIdentical(start(gr12), c(start(gr1), start(gr2)))
    checkIdentical(end(gr12), c(end(gr1), end(gr2)))
    checkIdentical(strand(gr12), c(strand(gr1), strand(gr2)))
    checkIdentical(mcols(gr12), rbind(mcols(gr1), mcols(gr2)))

    #########################################################################
    ## Concatenate GRanges objects with differing seqlevels
    x <- GRanges(seqinfo=Seqinfo(paste0("chr", 1:5), 1001:1005))
    y <- GRanges("chr4:1-11")
    current <- c(x, y)
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical("chr4", as.character(seqnames(current)))
    checkIdentical(seqlevels(x), levels(seqnames(current)))
    current <- c(y, x)
    checkTrue(validObject(current, complete=TRUE))
    checkIdentical("chr4", as.character(seqnames(current)))
    checkIdentical(paste0("chr", c(4, 1:3, 5)), levels(seqnames(current)))


    #########################################################################
    ## Concatenate GRanges objects with differing metadata columns
    colnames(mcols(gr2))[1] <- "score2"
    target <- c(gr1, gr2, ignore.mcols=TRUE)
    checkIdentical(dim(mcols(target)), c(length(target), 0L))
    mcols(target) <- rbind(cbind(mcols(gr1), score2=NA),
                           cbind(score=NA, mcols(gr2)))
    checkIdentical(target, c(gr1, gr2))

    #########################################################################
    ## More testing
    gr1 <- .make_TARGET_GRanges()
    gr2 <- gr1[, -1]
    target <- c(gr1, gr2, ignore.mcols=TRUE)
    mcols(target) <- rbind(mcols(gr1), cbind(score=NA, mcols(gr2)))
    checkIdentical(target, c(gr1, gr2))

    gr2 <- gr1
    names(gr2) <- NULL
    checkIdentical(as.data.frame(c(gr1, gr2), row.names=NULL),
                   rbind(as.data.frame(gr1, row.names=NULL),
                         as.data.frame(gr2, row.names=NULL)))
    checkIdentical(as.data.frame(c(gr2, gr1), row.names=NULL),
                   rbind(as.data.frame(gr2, row.names=NULL),
                         as.data.frame(gr1, row.names=NULL)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Split.
###

test_GRanges_split <- function()
{
    gr <- .make_TARGET_GRanges()
    checkException(split(gr, NULL), silent = TRUE)
    checkIdentical(split(gr, rep(c("a", "b"), each=5)),
                   GRangesList(a=head(gr, 5), b=tail(gr, 5)))
}

