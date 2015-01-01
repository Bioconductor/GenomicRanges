library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
cdsbytx <- cdsBy(txdb, "tx")[1:3]


### mapToTranscript and pmapToTranscript

test_mapToTranscript <- function()
{
    ## empty
    ans1 <- mapToTranscript(GRanges(), GRanges())
    ans2 <- mapToTranscript(GRanges(), GRangesList())
    checkTrue(length(ans1) == 0L)
    checkTrue(length(ans2) == 0L)
    checkIdentical(names(mcols(ans1)), c("xHits", "alignmentHits"))
    checkIdentical(names(mcols(ans2)), c("xHits", "alignmentHits"))

    ## strand
    x <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), 
                 width=200, names=LETTERS[1:3]), strand="+")
    align <- cdsbytx

    ans <- mapToTranscript(x, cdsbytx, ignore.strand=FALSE)
    checkTrue(length(ans) == 3L)
    checkIdentical(names(ans), c("B", "B", "C"))
    checkTrue(all(width(ans) == 200L))
    checkIdentical(mcols(ans)$xHits, c(2L, 2L, 3L))
    checkIdentical(mcols(ans)$alignmentHits, c(1L, 3L, 2L))

    strand(align[[2]][1]) <- "-"
    checkException(mapToTranscript(x, align, ignore.strand=FALSE), silent=TRUE) 

    strand(align[[2]]) <- "-"
    ans <- mapToTranscript(x, align, ignore.strand=FALSE) 
    checkIdentical(mcols(ans)$xHits, c(2L, 2L))
    checkIdentical(mcols(ans)$alignmentHits, c(1L, 3L))

    ## TODO: seqnames
}

test_mapToTranscript_range_order <- function()
{
    x <- GRanges("chrA", IRanges(43522349, width=1), strand="+")
    align1 <- GRangesList(GRanges("chrA", 
        IRanges(c(43522244, 43528406),
                c(43524145, 43528644)), strand="+"))
    align2 <- GRangesList(GRanges("chrA", 
        IRanges(c(43528406, 43522244),
                c(43528644, 43524145)), strand="+"))

    ## "+" strand
    ## smallest range first 
    ans <- mapToTranscript(x, align1, ignore.strand=FALSE)
    checkTrue(start(ans) == 106L)
    ## largest range first 
    ans <- mapToTranscript(x, align2, ignore.strand=FALSE)
    checkTrue(start(ans) == 106L)

    ## "-" strand
    strand(x) <- "-"
    strand(align1) <- "-"
    strand(align2) <- "-"
    ## smallest range first
    ans <- mapToTranscript(x, align1, ignore.strand=FALSE)
    checkTrue(start(ans) == 2036L)
    ## largest range first
    ans <- mapToTranscript(x, align2, ignore.strand=FALSE)
    checkTrue(start(ans) == 2036L)
}

test_pmapToTranscript <- function()
{
    ## empty
    ans1 <- pmapToTranscript(GRanges(), GRanges())
    ans2 <- pmapToTranscript(GRanges(), GRangesList())
    checkTrue(length(ans1) == 0L)
    checkTrue(length(ans2) == 0L)
    checkIdentical(names(mcols(ans1)), character(0))
    checkIdentical(names(mcols(ans2)), character(0))

    ## length
    x <- GRanges("chr1", IRanges(1, width=1))
    y <- GRanges("chr1", IRanges(1:5, width=1))
    checkException(pmapToTranscript(x, y), silent=TRUE)
    checkException(pmapToTranscript(x, GRangesList(y, y)), silent=TRUE)

    ## strand
    x <- GRanges("chr1", IRanges(c(6, 16, 1), width=1, names=LETTERS[1:3]))
    align <- GRangesList(GRanges("chr1", IRanges(5, 10)),
                         GRanges("chr1", IRanges(c(5, 15), width=6)),
                         GRanges("chr1", IRanges(5, 10)))

    strand(x) <- "-"
    strand(align) <- "+"
    ans <- pmapToTranscript(x, align, ignore.strand=FALSE)
    checkTrue(length(x) == length(x))
    checkTrue(all(width(ans) == 0L))
    checkIdentical(names(ans), LETTERS[1:3])

    strand(align[[2]][1]) <- "-"
    checkException(pmapToTranscript(x, align, ignore.strand=FALSE), silent=TRUE) 

    strand(align) <- "+"
    strand(x) <- "+"
    ans <- pmapToTranscript(x, align, ignore.strand=FALSE) 
    checkIdentical(width(ans), c(1L, 1L, 0L))
    checkIdentical(start(ans), c(2L, 8L, 1L))

    strand(align) <- "-"
    strand(x) <- "-"
    ans <- pmapToTranscript(x, align, ignore.strand=FALSE) 
    checkIdentical(width(ans), c(1L, 1L, 0L))
    checkIdentical(start(ans), c(5L, 5L, 1L))

    ## out of bounds
    x <- GRanges("chr1", IRanges(rep(40, 3), width=10))
    align <- GRanges("chr1", IRanges(c(1, 50, 100), width=20))
    ans <- pmapToTranscript(x, align) 
    checkIdentical(seqnames(ans), Rle(as.factor("unmapped"), 3)) 
}

### mapToGenome and pmapToGenome

test_mapToGenome <- function()
{
    ## empty
    ans1 <- mapToGenome(GRanges(), GRanges())
    ans2 <- mapToGenome(GRanges(), GRangesList())
    checkTrue(length(ans1) == 0L)
    checkTrue(length(ans2) == 0L)
    checkIdentical(names(mcols(ans1)), c("xHits", "alignmentHits"))
    checkIdentical(names(mcols(ans2)), c("xHits", "alignmentHits"))

    ## strand
    x <- GRanges("chr1", IRanges(rep(10, 5), width=2), strand="+")
    names(x) <- c("A", "B", "C", "foo", "A")
    gr <- GRanges("chr1", IRanges(c(10, 100), width=50), strand="+")
    align = GRangesList(gr, gr, gr, gr)
    names(align) <- c("C", "bar", "C", "A") 

    ans <- mapToGenome(x, align)
    checkTrue(length(ans) == 4L)
    checkIdentical(names(ans), c("A", "C", "C", "A"))
    checkTrue(all(width(ans) == 2L))
    checkIdentical(mcols(ans)$xHits, c(1L, 3L, 3L, 5L))
    checkIdentical(mcols(ans)$alignmentHits, c(4L, 1L, 3L, 4L))
    ans <- mapToGenome(x, unlist(align, use.names=TRUE))
    checkTrue(length(ans) == 8L)
    checkIdentical(names(ans), c("A", "A", rep("C", 4), "A", "A"))
    checkIdentical(mcols(ans)$alignmentHits, c(7L, 8L, 1L, 2L, 5L, 6L, 7L, 8L))

    strand(align[[1]][1]) <- "-"
    checkException(mapToGenome(x, align, ignore.strand=FALSE), silent=TRUE) 

    strand(align[[1]]) <- "-"
    ans <- mapToGenome(x, align, ignore.strand=FALSE) 
    checkIdentical(mcols(ans)$alignmentHits, c(4L, 3L, 4L))
    ans <- mapToGenome(x, align, ignore.strand=TRUE) 
    checkIdentical(mcols(ans)$alignmentHits, c(4L, 1L, 3L, 4L))

    ## TODO: seqnames

}

test_pmapToGenome <- function()
{
    ## empty
    ans1 <- pmapToGenome(GRanges(), GRanges())
    ans2 <- pmapToGenome(GRanges(), GRangesList())
    checkTrue(length(ans1) == 0L)
    checkTrue(length(ans2) == 0L)
    checkIdentical(names(mcols(ans1)), character(0))
    checkIdentical(names(mcols(ans2)), character(0))

    ## length
    x <- GRanges("chr1", IRanges(1, width=1))
    align <- GRanges("chr1", IRanges(1:5, width=1))
    checkException(pmapToGenome(x, align), silent=TRUE)
    checkException(pmapToGenome(x, GRangesList(align, align)), silent=TRUE)

    ## strand
    x <- GRanges("chr1", IRanges(c(1, 50, 150), width=1, names=LETTERS[1:3]))
    gr <- GRanges("chr1", IRanges(c(100, 300), width=100))
    align <- GRangesList(gr, gr, gr)

    strand(x) <- "-"
    strand(align) <- "+"
    ans <- pmapToGenome(x, align, ignore.strand=FALSE)
    checkTrue(length(x) == length(x))
    checkTrue(all(width(ans) == 0L))
    checkIdentical(names(ans), LETTERS[1:3])

    strand(align[[2]][1]) <- "-"
    checkException(pmapToGenome(x, align, ignore.strand=FALSE), silent=TRUE) 

    strand(align) <- "+"
    strand(x) <- "+"
    ans <- pmapToGenome(x, align, ignore.strand=FALSE) 
    checkIdentical(start(ans), c(100L, 149L, 349L))

    strand(align) <- "-"
    strand(x) <- "-"
    ans <- pmapToGenome(x, align, ignore.strand=FALSE) 
    checkIdentical(start(ans), c(399L, 350L, 150L))

    ## out of bounds
    x <- GRanges("chr1", IRanges(c(1, 50, 100), width=50))
    align <- GRanges("chr1", IRanges(rep(40, 3), width=20))
    ans <- pmapToGenome(x, align) 
    checkIdentical(seqnames(ans), Rle(as.factor("unmapped"), 3)) 
}
