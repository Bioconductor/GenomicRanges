group_id <- c("A", "B", "C", "C", "D", "D", "E", "F", "G", "H", "H")

features <- GRanges(
    seqnames = Rle(c("chr1", "chr2", "chr1", "chr1", "chr2", "chr2", 
        "chr1", "chr1", "chr2", "chr1", "chr1")),
    strand = strand(rep("+", length(group_id))),
    ranges = IRanges(
        start=c(1000, 2000, 3000, 3600, 7000, 7500, 4000, 4000, 3000, 
                5000, 5400),
    width=c(500, 900, 500, 300, 600, 300, 500, 900, 500, 500, 500)),
    DataFrame(group_id))

featuresList <- split(features, values(features)[["group_id"]])


reads <- GappedAlignments(
    names = c("a","b","c","d","e","f","g"),
    rname = Rle(c(rep(c("chr1", "chr2"), 3), "chr1")),
    pos = as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
    cigar = c("500M", "100M", "300M", "500M", "300M", 
        "50M200N50M", "50M150N50M"),
    strand = strand(rep.int("+", 7L)))

.getCounts <- function(res)
{
    as.integer(assays(res)$counts)
}


test_GRfeatures <- function()
{
    Union <- summarizeOverlaps(features, reads)
    IS <- summarizeOverlaps(features, reads, mode="IntersectionStrict")
    INE <- summarizeOverlaps(features, reads, mode="IntersectionNotEmpty")
 
    checkIdentical(.getCounts(Union), as.integer(c(1,1,0,0,0,0,0,0,1,0,0))) 
    checkIdentical(.getCounts(IS), as.integer(c(0,1,0,0,1,0,0,0,1,1,0))) 
    checkIdentical(.getCounts(INE), as.integer(c(1,1,0,0,1,0,0,0,1,1,0))) 
}


test_GRLfeatures <- function()
{
    Union <- summarizeOverlaps(featuresList, reads)
    IS <- summarizeOverlaps(featuresList, reads, mode="IntersectionStrict")
    INE <- summarizeOverlaps(featuresList, reads, mode="IntersectionNotEmpty")
 
    checkIdentical(.getCounts(Union), as.integer(c(1,1,1,1,0,0,1,1))) 
    checkIdentical(.getCounts(IS), as.integer(c(0,1,0,1,0,0,1,1))) 
    checkIdentical(.getCounts(INE), as.integer(c(1,1,1,1,0,0,1,1))) 
}


