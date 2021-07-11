make_test_GRangesList <- function() {
    a <- GRanges(Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
                 IRanges(1:10, end=10, names=head(letters, 10)),
                 strand=Rle(c("-", "+", "*", "+", "-"), c(1, 2, 2, 3, 2)),
                 seqinfo=Seqinfo(paste0("chr", 1:3)),
                 score=1:10, GC=seq(1, 0, length=10))
    b <- GRanges(Rle(c("chr2", "chr4", "chr5"), c(3, 6, 4)),
                 IRanges(1:13, end=13, names=tail(letters, 13)),
                 strand=Rle(c("-", "+", "-"), c(4, 5, 4)),
                 seqinfo=Seqinfo(paste0("chr", c(2, 4:5))),
                 score=1:13, GC=seq(0, 1, length=13))
    GRangesList(a=a, b=b)
}

test_GRangesList_construction <- function() {
    checkException(GRangesList(IRangesList()), silent = TRUE)

    checkTrue(validObject(new("CompressedGRangesList")))
    checkTrue(validObject(GRangesList()))
    checkTrue(validObject(GRangesList(GRanges())))
    checkTrue(validObject(GRangesList(GRanges(), GRanges())))
    checkTrue(validObject(GRangesList(list(GRanges(), GRanges()))))
    checkTrue(validObject(GRangesList(a = GRanges())))
    checkTrue(validObject(make_test_GRangesList()))
}

test_GRangesList_getters <- function() {
    grl <- make_test_GRangesList()
    checkIdentical(seqnames(grl), RleList(lapply(grl, seqnames), compress=TRUE))
    checkIdentical(ranges(grl), IRangesList(lapply(grl, ranges)))
    checkIdentical(strand(grl), RleList(lapply(grl, strand), compress=TRUE))
    checkIdentical(seqlengths(grl), seqlengths(grl@unlistData))
    checkIdentical(mcols(grl, level="within"),
                   SplitDataFrameList(lapply(grl, mcols)))
}

test_GRangesList_setters <- function() {
    grl0 <- GRangesList(A=GRanges("chr2", IRanges(3:2, 5)),
                        B=GRanges(c("chr2", "chrMT"), IRanges(7:6, 15)),
                        C=GRanges(c("chrY", "chrMT"), IRanges(17:16, 25)),
                        D=GRanges())

    current <- grl0
    seqlevels(current, pruning.mode="coarse") <- c("chr2", "chr5")
    target <- GRangesList(A=GRanges("chr2", IRanges(3:2, 5),
                                    seqinfo=Seqinfo(c("chr2", "chr5"))),
                          D=GRanges())
    checkIdentical(target, current)

    current <- grl0
    seqlevels(current, pruning.mode="fine") <- c("chr2", "chr5")
    target <- GRangesList(A=GRanges("chr2", IRanges(3:2, 5),
                                    seqinfo=Seqinfo(c("chr2", "chr5"))),
                          B=GRanges("chr2", IRanges(7, 15)),
                          C=GRanges(),
                          D=GRanges())
    checkIdentical(target, current)

    current <- grl0
    seqlevels(current, pruning.mode="tidy") <- c("chr2", "chr5")
    target <- GRangesList(A=GRanges("chr2", IRanges(3:2, 5),
                                    seqinfo=Seqinfo(c("chr2", "chr5"))),
                          B=GRanges("chr2", IRanges(7, 15)),
                          D=GRanges())
    checkIdentical(target, current)
}

test_GRangesList_coercion <- function() {
    ## as.data.frame
    gr1 <- GRanges(seqnames = c(1,1,2),
                   ranges = IRanges(1:3,4:6, names = head(letters,3)),
                   strand = strand(c("+", "-", "*")),
                   score = c(10L,2L,NA))
    gr2 <- GRanges(seqnames = c("chr1", "chr2"),
                   ranges = IRanges(1:2,1:2, names = tail(letters,2)),
                   strand = strand(c("*", "*")),
                   score = 12:13)
    grl <- GRangesList(a=gr1, b=gr2)
    df <-
      data.frame(group = togroup(PartitioningByWidth(grl)),
                 group_name = rep(c("a","b"), c(3, 2)),
                 seqnames = factor(c(1,1,2,"chr1","chr2")),
                 start = c(1:3,1:2), end = c(4:6,1:2),
                 width = c(4L, 4L, 4L, 1L, 1L),
                 strand = strand(c("+", "-", "*", "*", "*")),
                 score = c(10L,2L,NA,12:13),
                 stringsAsFactors = FALSE)
    checkIdentical(as.data.frame(grl), df)
}

test_GRangesList_IntegerRangesList <- function() {
    grl <- make_test_GRangesList()
    checkIdentical(start(grl), IntegerList(lapply(grl, start)))
    checkIdentical(end(grl), IntegerList(lapply(grl, end)))
    checkIdentical(width(grl), IntegerList(lapply(grl, width)))

    ## start
    checkException(start(grl) <- NULL, silent = TRUE)
    checkException(start(grl) <- 1:26, silent = TRUE)

    grl <- make_test_GRangesList()
    orig <- start(grl)
    start(grl) <- orig + 1L
    checkIdentical(start(grl), orig + 1L)

    ## end
    checkException(end(grl) <- NULL, silent = TRUE)
    checkException(end(grl) <- 1:26, silent = TRUE)

    grl <- make_test_GRangesList()
    orig <- end(grl)
    end(grl) <- orig + 1L
    checkIdentical(end(grl), orig + 1L)

    ## width
    checkException(width(grl) <- NULL, silent = TRUE)
    checkException(width(grl) <- 1:26, silent = TRUE)

    grl <- make_test_GRangesList()
    orig <- width(grl)
    width(grl) <- orig + 1L
    checkIdentical(width(grl), orig + 1L)
}

test_GRangesList_Vector <- function() {
    grl <- make_test_GRangesList()
    checkIdentical(grl, grl[])
    checkIdentical(grl[,"score"],
                   GRangesList(lapply(grl, function(x) x[,"score"])))
    checkIdentical(grl[seqnames(grl) == "chr2",],
                   GRangesList(lapply(grl, function(x)
                                      x[seqnames(x) == "chr2",])))
    checkIdentical(grl[seqnames(grl) == "chr2", "score"],
                   GRangesList(lapply(grl, function(x)
                                      x[seqnames(x) == "chr2", "score"])))

    checkIdentical(GRangesList(), c(GRangesList(), GRangesList()))
    checkIdentical(grl, c(grl, GRangesList()))
    checkIdentical(grl, c(GRangesList(), grl))
    GRL <- local({
        x <- grl
        names(x) <- toupper(names(grl))
        for (i in seq_len(length(x)))
            names(x[[i]]) <- toupper(names(grl[[i]]))
        x
    })
    res <- c(grl, GRL)
    checkTrue(validObject(res))
    ## [
    checkIdentical(grl, grl[Rle(TRUE)])
    checkIdentical(grl, res[seq_len(length(grl))])
    checkIdentical(GRL, res[-seq_len(length(grl))])
    ## checkException(c(grl, grl), "c() check for duplicated names", TRUE)

    checkIdentical(grl, local({
        x <- grl
        x[IRanges(1, 1)] <- grl[IRanges(1, 1)]
        x
    }), "[ by IRanges")
    checkIdentical(grl, local({
        x <- grl; x[1] <- grl[1]; x
    }), "[<- 1")
    checkIdentical(grl, local({
        x <- grl; x[1:2] <- grl[1:2]; x
    }), "[<- 2")
    checkIdentical(grl, local({
        x <- grl; x[[1]] <- grl[[1]]
        names(x) <- names(grl)
        x
    }), "[[<-, replace")
    checkIdentical(grl, local({
        x <- grl; x[["b"]] <- grl[[2]]
        names(x) <- names(grl)
        x
    }), "[[<-, replace char")
    checkIdentical(grl, local({
        x <- grl[1]; x[[2]] <- grl[[2]]
        names(x) <- names(grl)
        x
    }), "[[<-, extend-by-1")
    checkIdentical(grl, local({
        x <- grl[1]; x[["b"]] <- grl[[2]]
        names(x) <- names(grl)
        x
    }), "[[<-, extend-by-1 char")
}

