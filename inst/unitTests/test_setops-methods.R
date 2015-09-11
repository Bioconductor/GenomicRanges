###

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
    ## pintersect,GRanges,GRanges

    checkIdentical(GRanges(), pintersect(GRanges(), GRanges()))

    x <- GRanges("chr1", IRanges(c( 4, 13, 2,  9, 16, 9), width=1,
                                 names=letters[1:6]))
    checkIdentical(x, pintersect(x, x))
    checkIdentical(x, pintersect(x, unname(x)))
    checkIdentical(unname(x), pintersect(unname(x), x))

    current <- suppressWarnings(pintersect(x, GRanges("chrX", ranges(x))))
    checkTrue(all(width(current) == 0L))
    end(current) <- end(x)
    checkIdentical(x, current)

    y <- x
    strand(y) <- Rle(strand(c("*", "-", "+")), c(2, 2, 2))
    checkIdentical(y, pintersect(y, y))
    checkIdentical(y, pintersect(y, x))
    checkIdentical(y, pintersect(x, y))
    
    current <- pintersect(x, y, strict.strand=TRUE)
    same_strand_idx <- strand(x) == strand(y)
    checkTrue(all(width(current[!same_strand_idx]) == 0L))
    checkIdentical(x[same_strand_idx], current[same_strand_idx])
    end(current) <- end(x)
    checkIdentical(x, current)

    current <- pintersect(y, x, strict.strand=TRUE)
    same_strand_idx <- strand(x) == strand(y)
    checkTrue(all(width(current[!same_strand_idx]) == 0L))
    checkIdentical(y[same_strand_idx], current[same_strand_idx])
    end(current) <- end(y)
    checkIdentical(y, current)

    strand(x) <- rep(Rle(strand(c("*", "+", "-"))), 2)
    current <- pintersect(x, y)
    compat_strand_idx <- strand(x) == strand(y) |
                         strand(x) == "*" | strand(y) == "*"
    checkTrue(all(width(current[!compat_strand_idx]) == 0L))
    checkIdentical(ranges(x)[compat_strand_idx],
                   ranges(current)[compat_strand_idx])
    end(current) <- end(x)
    strand(x)[strand(x) == "*"] <- strand(y)[strand(x) == "*"]
    checkIdentical(x, current)

    checkIdentical(x, pintersect(x, y, ignore.strand=TRUE))
    checkIdentical(y, pintersect(y, x, ignore.strand=TRUE))

    x <- GRanges("chr1", IRanges(
             c( 2,  7,  7,  4, 13, 2,  9, 12,  9,  4,  8,  7,  5, 20),
             c( 6,  8, 12, 10, 12, 2, 12, 18,  8, 12, 11, 11, 11, 20)))

    y <- GRanges("chr1", IRanges(1, 20))
    checkIdentical(x, pintersect(x, y))

    y <- GRanges("chr1", IRanges(7, 11), strand="+")
    current <- pintersect(x, y)
    target <- x
    ranges(target) <- IRanges(
             c( 7,  7,  7,  7, 13, 2,  9, 12,  9,  7,  8,  7,  7, 20),
             c( 6,  8, 11, 10, 12, 1, 11, 11,  8, 11, 11, 11, 11, 19))
    nohit_idx <- c(5, 6, 14)
    strand(target)[-nohit_idx] <- strand(y)[-nohit_idx]
    checkIdentical(target, current)
 
    ## pintersect,GRangesList,GRanges and pintersect,GRanges,GRangesList

    checkIdentical(as(target, "GRangesList"),
                   pintersect(as(x, "GRangesList"), y))
    names(x) <- names(target) <- letters[1:14]
    checkIdentical(as(target, "GRangesList"),
                   pintersect(as(x, "GRangesList"), y))

    x <- GRangesList(x, rev(x))

    current <- pintersect(x, y)
    target <- mendoapply(pintersect, x, as(y, "GRangesList"))
    checkIdentical(target, current)

    y <- GRanges("chr1", IRanges(c(1, 7), c(20, 11)))
    strand(y) <- c("+", "-")

    current <- pintersect(x, y)
    target <- mendoapply(pintersect, x, as(y, "GRangesList"))
    checkIdentical(target, current)
    checkIdentical(target, pintersect(y, x))

    current <- reduce(pintersect(x, y, strict.strand=TRUE),
                      drop.empty.ranges=TRUE)
    target <- mendoapply(intersect, x, as(y, "GRangesList"))
    checkIdentical(target, current)

    ## pintersect,GRangesList,GRangesList

    gr <- make_test_GRanges()
    grlist <- make_test_GRangesList()
    expect <- reduce(grlist) 
    mcols(expect) <- NULL
    checkIdentical(pintersect(grlist, grlist), expect)
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
