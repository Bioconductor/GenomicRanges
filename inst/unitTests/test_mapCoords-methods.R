library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
gr <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), width = 200,
              names = LETTERS[1:3]))
cdsbytx <- cdsBy(txdb, "tx")

test_mapCoords <- function()
{
    x <- mapCoords(gr, cdsbytx)
    checkTrue(length(x) == 3L)
    checkIdentical(names(x), c("B", "B", "C"))
    checkIdentical(names(mcols(x)), c("queryHits", "subjectHits"))
    checkIdentical(start(x), c(645L, 609L, 1167L))
    checkTrue(all(width(x) == 200L)) 
}

test_mapCoords_strand <- function()
{
    gr <- GRanges("chrA", IRanges(c(43522349, 43522349),
                   width=c(1, 1)), strand=c("+", "-"))
    grl <- GRangesList(
               GRanges("chrA", IRanges(43528406, 43528644), strand="-"),
               GRanges("chrA", IRanges(43522244, 43524145), strand="-"))
    map <- mapCoords(gr, grl, ignore.strand=FALSE)
    checkIdentical(start(map) == 1797L)
    checkIdentical(end(map) == 1797L)

    map <- mapCoords(gr, grl, ignore.strand=TRUE)
    checkIdentical(start(map) == c(1797L, 1797L))
    checkIdentical(end(map) == c(1797L, 1797L))
}
