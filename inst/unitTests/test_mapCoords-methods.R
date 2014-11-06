library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
cdsbytx <- cdsBy(txdb, "tx")[1:3]

## mapCoords:

test_mapCoords_output <- function()
{
    cds <- cdsbytx
    from <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), 
                    width = 200, names = LETTERS[1:3]))
    x <- mapCoords(from, cds, ignore.strand=FALSE)
    checkTrue(length(x) == 3L)
    checkIdentical(names(x), c("B", "B", "C"))
    checkIdentical(names(mcols(x)), c("queryHits", "subjectHits"))
    checkTrue(all(width(x) == 200L))

    x <- mapCoords(from, cds, elt.hits=TRUE, ignore.strand=FALSE)
    checkTrue("eltHits" %in% names(mcols(x)))

    x <- mapCoords(from, cds, elt.loc=TRUE, ignore.strand=FALSE)
    checkTrue("eltLoc" %in% names(mcols(x)))
}

test_mapCoords_eltLoc <- function()
{
    cds <- cdsbytx
    from <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), width = 200)) 
    strand(from) <- "+" 
    x <- mapCoords(from, cds, elt.loc=TRUE, ignore.strand=FALSE)
    checkIdentical(start(x), c(645L, 609L, 1167L))
    checkIdentical(start(mcols(x)$eltLoc), c(208L, 172L, 333L))

    from <- GRanges("chr2L", IRanges(c(8200, 9000), width = 200))
    strand(from) <- "-" 
    strand(cds) <- "-"
    x <- mapCoords(from, cds, elt.loc=TRUE, ignore.strand=FALSE)
    checkIdentical(start(x), c(212L, 800L, 78L))
    checkIdentical(start(mcols(x)$eltLoc), c(212L, 191L, 78L))
}

test_mapCoords_range_order_pos <- function()
{
    from <- GRanges("chrA", IRanges(43522349, width=1), strand="+")
    ## Strand '+' smallest range first 
    grl <- GRangesList(GRanges("chrA", 
        IRanges(c(43522244, 43528406),
                c(43524145, 43528644)), strand="+"))
    x <- mapCoords(from, grl, ignore.strand=FALSE)
    checkTrue(start(x) == 106L)

    ## Strand '+' largest range first 
    grl <- GRangesList(GRanges("chrA", 
        IRanges(c(43528406, 43522244),
                c(43528644, 43524145)), strand="+"))
    x <- mapCoords(from, grl, ignore.strand=FALSE)
    checkTrue(start(x) == 106L)
}

test_mapCoords_range_order_neg <- function()
{
    from <- GRanges("chrA", IRanges(43522349, width=1), strand="-")

    ## Strand '-' smallest range first
    grl <- GRangesList(GRanges("chrA", 
        IRanges(c(43522244, 43528406),
                c(43524145, 43528644)), strand="-"))
    x <- mapCoords(from, grl, elt.loc = TRUE, ignore.strand=FALSE)
    checkTrue(start(x) == 2036L)
    checkTrue(start(mcols(x)$eltLoc) == 1797L)

    ## Strand '-' largest range first
    grl <- GRangesList(GRanges("chrA", 
        IRanges(c(43528406, 43522244),
                c(43528644, 43524145)), strand="-"))
    x <- mapCoords(from, grl, elt.loc = TRUE, ignore.strand=FALSE)
    checkTrue(start(x) == 2036L)
    checkTrue(start(mcols(x)$eltLoc) == 1797L)

    ## ignore.strand
    x <- mapCoords(from, grl, ignore.strand=TRUE)
    checkTrue(all(start(x) == 106L))
}

## pmapCoords:

test_pmapCoords <- function()
{
    cds <- cdsbytx
    from <- GRanges("chr2L", IRanges(c(7700, 8200, 8600, 8100), width = 200)) 
    checkException(pmapCoords(from, cds), silent=TRUE)

    x <- pmapCoords(from[1:3], cds, elt.loc=TRUE, elt.hits=TRUE)
    checkIdentical(start(x), c(21L, 445L))
    checkIdentical(start(mcols(x)$eltLoc), c(21L, 8L))
    checkIdentical(mcols(x)$queryHits, mcols(x)$subjectHits)

    gr1 <- GRanges("chr12", IRanges(c(30, 10), width=6), strand="+")
    gr2 <- GRanges("chr1", IRanges(c(20, 10, 100, 200), width=1), 
                   strand=c("-", "-", "+", "+"))
    to <- GRangesList(gr1, gr2)

    from <- GRanges("chr12", IRanges(c(30, 10), width=2))
    x <- pmapCoords(from, to, elt.loc=TRUE, elt.hits=TRUE)
    checkIdentical(start(x), 7L)
    checkIdentical(mcols(x)$eltHits, 1L)

    strand(from) <- "-"
    x <- pmapCoords(from, to, ignore.strand=FALSE)
    checkTrue(length(x) == 0L)

    from <- GRanges("chr1", IRanges(c(10, 100), width=1), strand=c("-", "+"))
    x <- pmapCoords(from, to, elt.loc=TRUE, elt.hits=TRUE, ignore.strand=FALSE)
    checkIdentical(mcols(x)$eltHits, 5L)
    from <- rev(from)
    x <- pmapCoords(from, to, elt.loc=TRUE, elt.hits=TRUE, ignore.strand=FALSE)
    checkIdentical(mcols(x)$eltHits, 4L)
}
