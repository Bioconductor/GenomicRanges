make_test_GRangesList <- function() {
    GRangesList(
        a = GRanges(
            seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
            ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
            strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
            seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")),
            score = 1:10, GC = seq(1, 0, length=10)),
        b = GRanges(
            seqnames = Rle(factor(c("chr2", "chr4", "chr5")), c(3, 6, 4)),
            ranges = IRanges(1:13, width = 13:1, names = tail(letters, 13)),
            strand = Rle(strand(c("-", "+", "-")), c(4, 5, 4)),
            seqinfo = Seqinfo(seqnames = paste("chr", c(2L, 4:5), sep="")),
            score = 1:13, GC = seq(0, 1, length=13))
    )
}

test_GRangesList_construction <- function() {
    checkException(GRangesList(IRangesList()), silent = TRUE)

    checkTrue(validObject(new("GRangesList")))
    checkTrue(validObject(GRangesList()))
    checkTrue(validObject(GRangesList(GRanges())))
    checkTrue(validObject(GRangesList(a = GRanges())))
    checkTrue(validObject(make_test_GRangesList()))
}

test_GRangesList_coercion <- function() {
    ## RangedDataList -> GRangesList
    rd <-
      RangedData(space = c(1,1,2),
                 ranges = IRanges(1:3,4:6, names = head(letters,3)),
                 strand = strand(c("+", "-", "*")),
                 score = c(10L,2L,NA))
    rdl <- RangedDataList(a = rd, b = rd)
    gr <-
      GRanges(seqnames = c(1,1,2),
              ranges = IRanges(1:3,4:6, names = head(letters,3)),
              strand = strand(c("+", "-", "*")),
              score = c(10L,2L,NA))
    grl <- GRangesList(a = gr, b = gr)
    checkIdentical(as(rdl, "GRangesList"), grl)

    ## as.data.frame
    gr1 <-
      GRanges(seqnames = c(1,1,2),
              ranges = IRanges(1:3,4:6, names = head(letters,3)),
              strand = strand(c("+", "-", "*")),
              score = c(10L,2L,NA))
    gr2 <-
      GRanges(seqnames = c("chr1", "chr2"),
              ranges = IRanges(1:2,1:2, names = tail(letters,2)),
              strand = strand(c("*", "*")),
              score = 12:13)
    grl <- GRangesList(a = gr1, b = gr2)
    df <-
      data.frame(group = togroup(grl),
                 group_name = rep(c("a","b"), c(3, 2)),
                 seqnames = factor(c(1,1,2,"chr1","chr2")),
                 start = c(1:3,1:2), end = c(4:6,1:2),
                 width = c(4L, 4L, 4L, 1L, 1L),
                 strand = strand(c("+", "-", "*", "*", "*")),
                 score = c(10L,2L,NA,12:13),
                 stringsAsFactors = FALSE)
    checkIdentical(as.data.frame(grl), df)
}

test_GRangesList_accessors <- function() {
    grl <- make_test_GRangesList()
    checkIdentical(seqnames(grl), RleList(lapply(grl, seqnames), compress=TRUE))
    checkIdentical(ranges(grl), IRangesList(lapply(grl, ranges)))
    checkIdentical(strand(grl), RleList(lapply(grl, strand), compress=TRUE))
    checkIdentical(seqlengths(grl), seqlengths(grl@unlistData))
    checkIdentical(mcols(grl, level="within"),
                   SplitDataFrameList(lapply(grl, mcols)))
}

test_GRangesList_RangesList <- function() {
    grl <- make_test_GRangesList()
    checkIdentical(start(grl), IntegerList(lapply(grl, start)))
    checkIdentical(end(grl), IntegerList(lapply(grl, end)))
    checkIdentical(width(grl), IntegerList(lapply(grl, width)))

    ## start
    checkException(start(GRangesList()) <- NULL, silent = TRUE)
    checkException(start(make_test_GRangesList()) <- 1:26, silent = TRUE)

    grl <- make_test_GRangesList()
    orig <- start(grl)
    start(grl) <- orig + 1L
    checkIdentical(start(grl), orig + 1L)

    ## end
    checkException(end(GRangesList()) <- NULL, silent = TRUE)
    checkException(end(make_test_GRangesList()) <- 1:26, silent = TRUE)

    grl <- make_test_GRangesList()
    orig <- end(grl)
    end(grl) <- orig + 1L
    checkIdentical(end(grl), orig + 1L)

    ## width
    checkException(width(GRangesList()) <- NULL, silent = TRUE)
    checkException(width(make_test_GRangesList()) <- 1:26, silent = TRUE)

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
