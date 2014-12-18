library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
cdsbytx <- cdsBy(txdb, "tx")[1:3]

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
    checkException(pmapToTranscript(x, GRanges()), silent=TRUE)
    checkException(pmapToTranscript(x, GRangesList()), silent=TRUE)

    ## strand
    x <- GRanges("chr1", IRanges(c(6, 16, 1), width=1))
    align <- GRangesList(GRanges("chr1", IRanges(5, 10)),
                         GRanges("chr1", IRanges(c(5, 15), width=6)),
                         GRanges("chr1", IRanges(5, 10)))

    strand(x) <- "-"
    strand(align) <- "+"
    ans <- pmapToTranscript(x, align, ignore.strand=FALSE)
    checkTrue(length(x) == length(x))
    checkTrue(all(width(ans) == 0L))

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

    ## TODO: seqnames
}
