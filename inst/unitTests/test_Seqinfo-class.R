gr1 <- GRanges(c("chr1", "chr1"), IRanges(1, 5))
gr2 <- GRanges(c("chr1", "chr2"), IRanges(1, 5))
gr3 <- GRanges(c("chr3", "chr3", "chr3"), IRanges(1:3, 5))

test_Seqinfo_seqlevels_subset <- function()
{
    ## GRanges
    gr <- gr2
    seqlevels(gr, force=TRUE) <- "chr2"
    checkIdentical("chr2", seqlevels(gr))
    checkIdentical(1L, length(gr))
    gr <- suppressWarnings(c(gr1, gr2, gr3))

    ## GRangesList
    grl <- GRangesList(gr1, gr2, gr3)
    seqlevels(grl, force=TRUE) <- "chr1"
    checkIdentical("chr1", seqlevels(grl))
    checkIdentical(1L, length(grl))

    grl <- GRangesList(gr2)
    seqlevels(grl, force=TRUE) <- "chr1"
    checkIdentical(0L, length(grl)) ## warning?

    ## GAlignments
    gal <- GAlignments(seqnames=Rle(c("chr1", "chr2")),
                        pos=as.integer(c(10, 100)),
                        cigar=c("50M", "50M"),
                        strand=strand(c("*", "*")))
    seqlevels(gal, force=TRUE) <- "chr2"
    checkIdentical("chr2", seqlevels(gal))
}

test_Seqinfo_seqlevels_rename <- function()
{
    ## GRanges
    gr <- c(gr1, gr2, gr3) 
    seqlevels(gr) <- gsub("chr", "CHR", seqlevels(gr))
    checkIdentical(c("CHR1", "CHR2", "CHR3"), seqlevels(gr))
    seqlevels(gr)[seqlevels(gr) == "CHR2"] <- "2"
    checkIdentical(c("CHR1", "2", "CHR3"), seqlevels(gr))
    #checkException(seqlevels(gr) <- rev(seqlevels(gr))) ## exception?
    #nms <- seqlevels(gr)
    #new <- rev(seqlevels(gr))
    #names(new) <- nms
    #checkException(seqlevels(gr) <- rev(seqlevels(gr))) ## exception?

    ## GRangesList
    grl <- GRangesList(gr1, gr2, gr3)
    #seqlevels(grl)[seqlevels(grl) == "chr2"] <- "chr1" ## should work?
    #seqlevels(grl)[seqlevels(grl) == "chr3"] <- "chr1" ## should work?
    seqlevels(grl)[seqlevels(grl) == "chr3"] <- "3"
    checkIdentical(c("chr1", "chr2", "3"), seqlevels(grl))


    ## GAlignments
    gal <- GAlignments(seqnames=Rle(c("chr1", "chr2")),
                        pos=as.integer(c(10, 100)),
                        cigar=c("50M", "50M"),
                        strand=strand(c("*", "*")))
    seqlevels(gal)[seqlevels(gal) == "chr2"] <- "2"
    checkIdentical(c("chr1", "2"),  seqlevels(gal))
}

test_Seqinfo.merge <- function()
{
    x <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chrM"),
                 seqlengths=c(100, 200, NA, 15),
                 isCircular=c(NA, FALSE, FALSE, TRUE))
    checkIdentical(merge(x), x)
    checkIdentical(merge(x, NULL), x)
    checkIdentical(merge(NULL, x), x)
    checkIdentical(merge(x, Seqinfo()), x)
    checkIdentical(merge(Seqinfo(), x), x)
    checkIdentical(merge(x, x), x)
    checkIdentical(merge(x, Seqinfo(rev(names(seqlengths)))), x)

    y <- Seqinfo(seqnames=c("chrM", "chr4", "chr3"),
                 seqlengths=c(15, NA, 300))
    got <- merge(x, y)
    want <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chrM", "chr4"),
                    seqlengths=c(100, 200, 300, 15, NA),
                    isCircular=c(NA, FALSE, FALSE, TRUE, NA))
    checkIdentical(got, want)

    ## This contradicts what 'x' says about circularity of chr3 and chrM:
    isCircular(y)[c("chr3", "chrM")] <- c(TRUE, FALSE)
    checkException(merge(x, y), silent=TRUE)
}

