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

    ## Strand '-' out of order
    grl <- GRangesList(GRanges("chrA", 
        IRanges(c(43522244, 43528406),
                c(43524145, 43528644)), strand="-"))
    map_neg_out <- mapCoords(gr, grl, ignore.strand=FALSE)
    checkIdentical(start(map_neg_out), start(map_neg_out))

    ## Strand '-' in order
    grl <- GRangesList(GRanges("chrA", 
        IRanges(c(43528406, 43522244),
                c(43528644, 43524145)), strand="-"))
    map_neg_in <- mapCoords(gr, grl, ignore.strand=FALSE)
    checkTrue(start(map_neg_in) == 1797L)
    checkTrue(end(map_neg_in) == 1797L)

    ## Strand '+' out of order
    grl <- GRangesList(GRanges("chrA", 
        IRanges(c(43528406, 43522244),
                c(43528644, 43524145)), strand="+"))
    map_pos_out <- mapCoords(gr, grl, ignore.strand=FALSE)
    checkTrue(start(map_pos_out) == 106L)

    ## Strand '+' in order
    grl <- GRangesList(GRanges("chrA", 
        IRanges(c(43522244, 43528406),
                c(43524145, 43528644)), strand="+"))
    map_pos_in <- mapCoords(gr, grl, ignore.strand=FALSE)
    checkTrue(start(map_pos_out) == 106L)

    ## ignore.strand
    map_pos_is <- mapCoords(gr, grl, ignore.strand=TRUE)
    checkTrue(all(start(map_pos_is) == 106L))
}
