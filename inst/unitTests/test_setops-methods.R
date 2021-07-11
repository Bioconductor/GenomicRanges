###

make_test_GRanges <- function(seqlevels=paste0("chr", 1:5))
    GRanges(Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
            IRanges(1:10, end=10, names=head(letters, 10)),
            Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
            seqinfo=Seqinfo(seqlevels),
            score=1:10, GC=seq(1, 0, length=10))

make_test_GRangesList <- function() {
    a <- make_test_GRanges(paste0("chr", 1:3))
    b <- GRanges(Rle(factor(c("chr2", "chr4", "chr5")), c(3, 6, 4)),
                 IRanges(1:13, end=13, names=tail(letters, 13)),
                 Rle(strand(c("-", "+", "-")), c(4, 5, 4)),
                 seqinfo=Seqinfo(paste0("chr", c(2, 4:5))),
                 score=1:13, GC=seq(0, 1, length=13))
    GRangesList(a=a, b=b)
}

test_union <- function()
{
    gr <- make_test_GRanges()
    grlist <- make_test_GRangesList()

    checkException(union(gr, grlist), silent = TRUE)
    checkException(union(grlist, gr), silent = TRUE)

    expect <- gr
    mcols(expect) <- NULL
    expect <- reduce(gr)
    checkIdentical(union(gr, gr), expect)

    checkIdentical(union(grlist[[1]], grlist[[2]]), reduce(unlist(grlist)))
}

test_intersect <- function()
{
    gr <- make_test_GRanges()
    grlist <- make_test_GRangesList()

    checkException(intersect(gr, grlist), silent = TRUE)
    checkException(intersect(grlist, gr), silent = TRUE)

    expect <- gr
    mcols(expect) <- NULL
    expect <- reduce(gr)
    checkIdentical(intersect(gr, gr), expect)

    checkIdentical(union(grlist[[1]], grlist[[2]]), reduce(unlist(grlist)))

    gr1 <-
      GRanges(seqnames = c("chr1", "chr2", "chr3", "chr4"),
              ranges = IRanges(1:4, width = 10),
              strand = c("-", "+", "+", "+"))
    gr2 <-
      GRanges(seqnames = c("chr1", "chr2", "chr5", "chr6"),
              ranges = IRanges(1:4, width = 10),
              strand = c("-", "+", "+", "+"))
    expect <-
      GRanges(seqnames =
              factor(c("chr1", "chr2"), levels = paste("chr", 1:6, sep = "")),
              ranges = IRanges(1:2, width = 10),
              strand = c("-", "+"))
    checkIdentical(suppressWarnings(intersect(gr1, gr2)), expect)

    ## intersect,GRangesList,GRangesList

    gr <- make_test_GRanges()
    grlist <- make_test_GRangesList()
    expect <- reduce(grlist)
    mcols(expect) <- NULL
    checkIdentical(intersect(grlist, grlist), expect)
}

test_setdiff <- function()
{
    gr <- make_test_GRanges()
    grlist <- make_test_GRangesList()

    checkException(setdiff(gr, grlist), silent = TRUE)
    checkException(setdiff(grlist, gr), silent = TRUE)

    expect <- GRanges(seqlengths=seqlengths(gr))
    checkIdentical(setdiff(gr, gr), expect)

    checkIdentical(setdiff(grlist[[1]], grlist[[2]]), reduce(grlist[[1]]))

    gr1 <-
      GRanges(seqnames = c("chr1", "chr2", "chr3", "chr4"),
              ranges = IRanges(1:4, width = 10),
              strand = c("-", "+", "+", "+"))
    gr2 <-
      GRanges(seqnames = c("chr1", "chr2", "chr5", "chr6"),
              ranges = IRanges(1:4, width = 10),
              strand = c("-", "+", "+", "+"))
    expect <-
      GRanges(seqnames =
              factor(c("chr3", "chr4"), levels = paste("chr", 1:6, sep = "")),
              ranges = IRanges(3:4, width = 10),
              strand = c("+", "+"))
    checkIdentical(suppressWarnings(setdiff(gr1, gr2)), expect)
}

test_punion <- function()
{
    gr <- make_test_GRanges()
    grlist <- make_test_GRangesList()

    expect <- gr
    mcols(expect) <- NULL
    checkIdentical(punion(gr, gr), expect)
}

test_pintersect <- function()
{
    ## pintersect,GRanges,GRanges

    gr0 <- GRanges()
    current <- pintersect(gr0, gr0, drop.nohit.ranges=TRUE)
    target <- gr0
    checkIdentical(target, current)

    current <- pintersect(gr0, gr0)
    mcols(target)$hit <- logical(0)
    checkIdentical(target, current)

    x <- GRanges("chr1", IRanges(c( 4, 13, 2,  9, 16, 9), width=1,
                                 names=letters[1:6]))
    current <- pintersect(x, x)
    target <- x
    mcols(target)$hit <- TRUE
    checkIdentical(target, current)
    checkIdentical(target, pintersect(x, unname(x)))
    checkIdentical(unname(target), pintersect(unname(x), x))

    current <- suppressWarnings(pintersect(x, GRanges("chrX", ranges(x))))
    target <- x
    mcols(target)$hit <- FALSE
    width(target) <- 0
    checkIdentical(target, current)

    current <- suppressWarnings(pintersect(x, GRanges("chrX", ranges(x)),
                                           drop.nohit.ranges=TRUE))
    checkIdentical(x[0], current)

    y <- x
    strand(y) <- Rle(strand(c("*", "-", "+")), c(2, 2, 2))
    target <- y
    mcols(target)$hit <- TRUE
    checkIdentical(target, pintersect(y, y))
    checkIdentical(target, pintersect(y, x))
    checkIdentical(target, pintersect(x, y))

    current <- pintersect(x, y, strict.strand=TRUE)
    target <- x
    mcols(target)$hit <- as.logical(strand(x) == strand(y))
    width(target)[!mcols(target)$hit] <- 0
    checkIdentical(target, current)

    current <- pintersect(y, x, strict.strand=TRUE)
    strand(target) <- strand(y)
    checkIdentical(target, current)

    strand(x) <- rep(Rle(strand(c("*", "+", "-"))), 2)
    current <- pintersect(x, y)
    target <- x
    mcols(target)$hit <- as.logical(strand(x) == strand(y) |
                                    strand(x) == "*" | strand(y) == "*")
    width(target)[!mcols(target)$hit] <- 0
    disambig_strand_idx <- strand(target) == "*" & mcols(target)$hit
    strand(target)[disambig_strand_idx] <- strand(y)[disambig_strand_idx]
    checkIdentical(target, current)

    current <- pintersect(x, y, ignore.strand=TRUE)
    target <- x
    mcols(target)$hit <- TRUE
    checkIdentical(target, current)

    current <- pintersect(y, x, ignore.strand=TRUE)
    target <- y
    mcols(target)$hit <- TRUE
    checkIdentical(target, current)

    x <- GRanges("chr1", IRanges(
             c( 2,  7,  7,  4, 13, 2,  9, 12,  9,  4,  8,  7,  5, 20),
             c( 6,  8, 12, 10, 12, 2, 12, 18,  8, 12, 11, 11, 11, 20)))

    y <- GRanges("chr1", IRanges(1, 20))
    current <- pintersect(x, y)
    target <- x
    mcols(target)$hit <- TRUE
    checkIdentical(target, current)

    y <- GRanges("chr1", IRanges(7, 11), strand="+")
    current <- pintersect(x, y)
    target <- x
    ranges(target) <- IRanges(
             c( 7,  7,  7,  7, 13, 2,  9, 12,  9,  7,  8,  7,  7, 20),
             c( 6,  8, 11, 10, 12, 1, 11, 11,  8, 11, 11, 11, 11, 19))
    mcols(target)$hit <- TRUE
    nohit_idx <- c(5, 6, 14)
    mcols(target)$hit <- !(1:14 %in% nohit_idx)
    strand(target)[mcols(target)$hit] <- strand(y)
    checkIdentical(target, current)

    current <- pintersect(x, y, drop.nohit.ranges=TRUE)
    target <- target[mcols(target)$hit]
    mcols(target)$hit <- NULL
    checkIdentical(target, current)

    ## pintersect,GRangesList,GRanges and pintersect,GRanges,GRangesList

    x <- GRangesList(x, rev(x))

    current <- pintersect(x, y)
    target <- mendoapply(pintersect, x, as(y, "CompressedGRangesList"))
    checkIdentical(target, current)

    current <- pintersect(x, y, drop.nohit.ranges=TRUE)
    target <- mendoapply(pintersect, x, as(y, "CompressedGRangesList"),
                         MoreArgs=list(drop.nohit.ranges=TRUE))
    checkIdentical(target, current)

    y <- GRanges("chr1", IRanges(c(1, 7), c(20, 11)))
    strand(y) <- c("+", "-")

    current <- pintersect(x, y)
    target <- mendoapply(pintersect, x, as(y, "CompressedGRangesList"))
    checkIdentical(target, current)
    checkIdentical(target, pintersect(y, x))

    current <- pintersect(x, y, drop.nohit.ranges=TRUE)
    target <- mendoapply(pintersect, x, as(y, "CompressedGRangesList"),
                         MoreArgs=list(drop.nohit.ranges=TRUE))
    checkIdentical(target, current)

    strand(x[[1]]) <- c("+", "-")

    current <- pintersect(x, y)
    target <- mendoapply(pintersect, x, as(y, "CompressedGRangesList"))
    checkIdentical(target, current)
    checkIdentical(target, pintersect(y, x))

    current <- pintersect(x, y, drop.nohit.ranges=TRUE)
    target <- mendoapply(pintersect, x, as(y, "CompressedGRangesList"),
                         MoreArgs=list(drop.nohit.ranges=TRUE))
    checkIdentical(target, current)
    checkIdentical(target, pintersect(y, x, drop.nohit.ranges=TRUE))

    current <- pintersect(x, y, strict.strand=TRUE)
    target <- mendoapply(pintersect, x, as(y, "CompressedGRangesList"),
                         MoreArgs=list(strict.strand=TRUE))
    checkIdentical(target, current)
    checkIdentical(target, pintersect(y, x, strict.strand=TRUE))

    current <- pintersect(x, y, drop.nohit.ranges=TRUE, strict.strand=TRUE)
    target <- mendoapply(pintersect, x, as(y, "CompressedGRangesList"),
                         MoreArgs=list(drop.nohit.ranges=TRUE,
                                       strict.strand=TRUE))
    checkIdentical(target, current)
    checkIdentical(target, pintersect(y, x, drop.nohit.ranges=TRUE,
                                            strict.strand=TRUE))

    current <- reduce(pintersect(x, y, strict.strand=TRUE),
                      drop.empty.ranges=TRUE)
    target <- mendoapply(intersect, x, as(y, "CompressedGRangesList"))
    checkIdentical(target, current)

    current <- reduce(pintersect(x, y, ignore.strand=TRUE),
                      drop.empty.ranges=TRUE, ignore.strand=TRUE)
    target <- mendoapply(intersect, x, as(y, "CompressedGRangesList"),
                         MoreArgs=list(ignore.strand=TRUE))
    checkIdentical(target, current)
}

test_psetdiff <- function()
{
    gr <- make_test_GRanges()
    grlist <- make_test_GRangesList()

    checkException(psetdiff(grlist, gr), silent = TRUE)

    expect <- gr
    width(expect) <- 0L
    mcols(expect) <- NULL
    checkIdentical(psetdiff(gr, gr), expect)

    expect <- gr
    mcols(expect) <- NULL
    strand(expect) <- Rle(c("-", "+", "-"), c(1, 7, 2))
    end(expect)[5:6] <- c(5, 5)
    checkIdentical(psetdiff(gr, rev(gr)), expect)

    expect <-
      GRangesList(a = GRanges(),
                  b =
                  GRanges(factor("chr2",
                                 levels = paste("chr", 1:5, sep = "")),
                          IRanges(2, 10),
                          "+"))
    checkIdentical(psetdiff(gr[1:2], grlist[1:2]), expect)
}
