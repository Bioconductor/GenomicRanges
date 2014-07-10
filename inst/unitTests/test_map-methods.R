
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
gr <- GRanges("chr2L", IRanges(c(7500, 8400, 9000), width = 200,
              names = LETTERS[1:3]))
cdsbytx <- cdsBy(txdb, "tx")

test_map_GRangesMapping <- function()
{
    x <- map(gr, cdsbytx)
    checkTrue(length(x) == 2L)
    checkTrue(length(hits(x)) == length(granges(x)))
    checkTrue(all(width(granges(x) == 200L))) 
}

test_map <- function()
{
    x <- map(gr, cdsbytx)
    checkTrue(all(names(granges(x)) %in% c("B", "C")))
    checkTrue(start(granges(x)) == c(645, 609, 1167))

    xx <- as(x, "GenomicRanges")
    checkTrue(all(names(mcols(xx)) %in% c("queryHits", "subjectHits"))
}




