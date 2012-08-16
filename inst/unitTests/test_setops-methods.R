make_test_GRanges <- function() {
    new("GRanges",
        seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3"),
                              levels = paste("chr", 1:5, sep = "")),
                       c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
        seqinfo = Seqinfo(seqnames = paste("chr", 1:5, sep="")),
        elementMetadata = DataFrame(score = 1:10, GC = seq(1, 0, length=10)))
}

make_test_GRangesList <- function() {
    suppressWarnings(GRangesList(
        a =
        new("GRanges",
            seqnames = Rle(factor(c("chr1", "chr2", "chr1", "chr3")), c(1, 3, 2, 4)),
            ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
            strand = Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
            seqinfo = Seqinfo(seqnames = paste("chr", 1:3, sep="")),
            elementMetadata = DataFrame(score = 1:10, GC = seq(1, 0, length=10))),
        b =
        new("GRanges",
            seqnames = Rle(factor(c("chr2", "chr4", "chr5")), c(3, 6, 4)),
            ranges = IRanges(1:13, width = 13:1, names = tail(letters, 13)),
            strand = Rle(strand(c("-", "+", "-")), c(4, 5, 4)),
            seqinfo = Seqinfo(seqnames = paste("chr", c(2L, 4:5), sep="")),
            elementMetadata = DataFrame(score = 1:13, GC = seq(0, 1, length=13)))))
}

test_union <- function()
{
    gr <- make_test_GRanges()
    grlist <- make_test_GRangesList()

    checkException(union(gr, grlist), silent = TRUE)
    checkException(union(grlist, gr), silent = TRUE)
    checkException(union(grlist, grlist), silent = TRUE)

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
    checkException(intersect(grlist, grlist), silent = TRUE)

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
}

test_setdiff <- function()
{
    gr <- make_test_GRanges()
    grlist <- make_test_GRangesList()

    checkException(setdiff(gr, grlist), silent = TRUE)
    checkException(setdiff(grlist, gr), silent = TRUE)
    checkException(setdiff(grlist, grlist), silent = TRUE)

    expect <- GRanges(seqnames(gr)[integer(0)], seqlengths = seqlengths(gr))
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

    checkException(punion(grlist, grlist), silent = TRUE)

    expect <- gr
    mcols(expect) <- NULL
    checkIdentical(punion(gr, gr), expect)
}

test_pintersect <- function()
{
    gr <- make_test_GRanges()
    grlist <- make_test_GRangesList()

    expect <- reduce(grlist) 
    mcols(expect) <- NULL
    checkIdentical(pintersect(grlist, grlist), expect)

    expect <- gr
    mcols(expect) <- NULL
    checkIdentical(pintersect(gr, gr), expect)
    checkIdentical(pintersect(gr[1:2], grlist), pintersect(grlist, gr[1:2]))

    expect <-
      GRangesList(a =
                  GRanges(factor(c("chr2", "chr2", "chr2"),
                                 levels = paste("chr", 1:5, sep = "")),
                          IRanges(c(4,4,4), c(10,10,10), names = letters[2:4]),
                          c("+", "+", "*")),
                  b = GRanges())
    checkIdentical(pintersect(gr[4:5], grlist), expect)
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
