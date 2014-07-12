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
