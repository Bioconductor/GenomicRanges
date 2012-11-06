.getCounts <- function(res)
{
    as.integer(assays(res)$counts)
}

test_summarizeOverlaps_Union_single <- function()
{
    ## single-end no gaps
    mode <- "Union"
    ga <- GappedAlignments("chr1", 20L, "11M", strand("+")) 
    ann <- GRanges("chr1", IRanges(c(1, 10, 25, 22), c(50, 25, 40, 26)), "+")
    res <- summarizeOverlaps(ann[1], ga, mode)
    checkIdentical(1L, .getCounts(res)) 
    res <- summarizeOverlaps(ann[2], ga, mode)
    checkIdentical(1L, .getCounts(res)) 
    res <- summarizeOverlaps(ann[3], ga, mode)
    checkIdentical(1L, .getCounts(res)) 
    res <- summarizeOverlaps(ann[4], ga, mode)
    checkIdentical(1L, .getCounts(res)) 
    ## >1 feature
    res <- summarizeOverlaps(ann, ga, mode)
    checkIdentical(c(0L, 0L, 0L, 0L), .getCounts(res)) 
}

test_summarizeOverlaps_Union_paired <- function()
{
    ## single-end with gaps behave like paired-end
    mode <- "Union"
    ga <- GappedAlignments("chr1", 1L, "10M4N11M", strand("+"))
    ga1 <- GappedAlignments("chr1", 1L, "10M", strand("+"))
    ga2 <- GappedAlignments("chr1", 15L, "11M", strand("-"))
    galp <- GappedAlignmentPairs(ga1, ga2, TRUE)
    ann <- GRanges("chr1", IRanges(c(1, 5, 12, 20), c(25, 20, 14, 30)), "+")

    res_ga <- summarizeOverlaps(ann[1], ga, mode)
    res_galp <- summarizeOverlaps(ann[1], galp, mode)
    checkIdentical(1L, .getCounts(res_ga))
    checkIdentical(1L, .getCounts(res_galp))
    res_ga <- summarizeOverlaps(ann[2], ga, mode)
    res_galp <- summarizeOverlaps(ann[2], galp, mode)
    checkIdentical(1L, .getCounts(res_ga))
    checkIdentical(1L, .getCounts(res_galp))
    res_ga <- summarizeOverlaps(ann[3], ga, mode)
    res_galp <- summarizeOverlaps(ann[3], galp, mode)
    checkIdentical(0L, .getCounts(res_ga))
    checkIdentical(0L, .getCounts(res_galp))
    res_ga <- summarizeOverlaps(ann[4], ga, mode)
    res_galp <- summarizeOverlaps(ann[4], galp, mode)
    checkIdentical(1L, .getCounts(res_ga))
    checkIdentical(1L, .getCounts(res_galp))
    ## >1 feature
    res_ga <- summarizeOverlaps(ann, ga, mode)
    res_galp <- summarizeOverlaps(ann, ga, mode)
    checkIdentical(c(0L, 0L, 0L, 0L), .getCounts(res_ga))
    checkIdentical(c(0L, 0L, 0L, 0L), .getCounts(res_galp))
}

test_summarizeOverlaps_IntersectionStrict_single <- function()
{
    ## single-end, no gaps 
    mode <- "IntersectionStrict"
    ga <- GappedAlignments("chr1", 7L, "6M", strand("+")) 
    ann <- GRanges("chr1", IRanges(c(1, 5, 10), width=10), "+") 
    res <- summarizeOverlaps(ann[1], ga, mode)
    checkIdentical(0L, .getCounts(res)) 
    res <- summarizeOverlaps(ann[2], ga, mode)
    checkIdentical(1L, .getCounts(res)) 
    res <- summarizeOverlaps(ann[3], ga, mode)
    checkIdentical(0L, .getCounts(res))
    ## >1 feature 
    ann <- GRanges("chr1", IRanges(c(5, 6, 10), c(15, 16, 15)), "+") 
    res <- summarizeOverlaps(ann[1:2], ga, mode)
    checkIdentical(c(0L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ann[c(1,3)], ga, mode)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
}

test_summarizeOverlaps_IntersectionStrict_paired <- function()
{
    ## single-end with gaps behave like paired-end
    mode <- "IntersectionStrict"
    ga <- GappedAlignments("chr1", 10L, "6M4N6M", strand("+"))
    ga1 <- GappedAlignments("chr1", 10L, "6M", strand("+"))
    ga2 <- GappedAlignments("chr1", 20L, "6M", strand("-"))
    galp <- GappedAlignmentPairs(ga1, ga2, TRUE)
    ann <- GRanges("chr1", IRanges(c(1, 1, 20), c(30, 15, 30)), "+")

    res_ga <- summarizeOverlaps(ann[1], ga, mode)
    res_galp <- summarizeOverlaps(ann[1], galp, mode)
    checkIdentical(1L, .getCounts(res_ga))
    checkIdentical(1L, .getCounts(res_galp)) 
    res_ga <- summarizeOverlaps(ann[2], ga, mode)
    res_galp <- summarizeOverlaps(ann[2], galp, mode)
    checkIdentical(0L, .getCounts(res_ga))
    checkIdentical(0L, .getCounts(res_galp)) 
    res_ga <- summarizeOverlaps(ann[3], ga, mode)
    res_galp <- summarizeOverlaps(ann[3], galp, mode)
    checkIdentical(0L, .getCounts(res_ga))
    checkIdentical(0L, .getCounts(res_galp)) 
    ## >1 feature
    res_ga <- summarizeOverlaps(ann, ga, mode)
    res_galp <- summarizeOverlaps(ann, galp, mode)
    checkIdentical(c(1L, 0L, 0L), .getCounts(res_ga))
    checkIdentical(c(1L, 0L, 0L), .getCounts(res_galp))
}

test_summarizeOverlaps_IntersectionNotEmpty_single <- function()
{
    ## single-end, no gaps
    mode <- "IntersectionNotEmpty"
    ga <- GappedAlignments("chr1", 10L, "11M", strand("+")) 
    ann <- GRanges("chr1", IRanges(c(1, 5, 12), c(15, 30, 15)), "+") 
    res <- summarizeOverlaps(ann[1], ga, mode)
    checkIdentical(1L, .getCounts(res)) 
    res <- summarizeOverlaps(ann[2], ga, mode)
    checkIdentical(1L, .getCounts(res)) 
    res <- summarizeOverlaps(ann[3], ga, mode)
    checkIdentical(1L, .getCounts(res)) 
    res <- summarizeOverlaps(ann, ga, mode)
    checkIdentical(c(0L, 1L, 0L), .getCounts(res))
    ## >1 feature 
    ann <- GRanges("chr1", IRanges(c(5, 15), c(20, 25)), "+") 
    res <- summarizeOverlaps(ann, ga, mode)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
    ann <- GRanges("chr1", IRanges(c(5, 12), c(18, 25)), "+") 
    res <- summarizeOverlaps(ann, ga, mode)
    checkIdentical(c(0L, 0L), .getCounts(res))

    ann <- GRanges(rep("chr1", 3), IRanges(c(1L, 20L, 20L), 
        width=c(50, 11, 11)), c("+", "+", "-")) 
    ga <- GappedAlignments("chr1", 23L, "5M", strand("*"))
    res <- summarizeOverlaps(ann, ga, mode)
    checkIdentical(c(0L, 0L, 0L), .getCounts(res))
    strand(ga) <- "-"
    res <- summarizeOverlaps(ann, ga, mode)
    checkIdentical(c(0L, 0L, 1L), .getCounts(res))
    ga <- GappedAlignments("chr1", 28L, "5M", strand("+"))
    res <- summarizeOverlaps(ann, ga, mode)
    checkIdentical(c(1L, 0L, 0L), .getCounts(res))
}

test_summarizeOverlaps_IntersectionNotEmpty_paired <- function()
{
    ## single-end with gaps behave like paired-end
    mode <- "IntersectionNotEmpty"
    ga <- GappedAlignments("chr1", 10L, "6M4N6M", strand("+"))
    ga1 <- GappedAlignments("chr1", 10L, "6M", strand("+"))
    ga2 <- GappedAlignments("chr1", 20L, "6M", strand("-"))
    galp <- GappedAlignmentPairs(ga1, ga2, TRUE)
    ann <- GRanges("chr1", IRanges(c(1, 1, 20), c(30, 15, 30)), "+")

    ## single-end, gaps 
    res_ga <- summarizeOverlaps(ann[1], ga, mode)
    res_galp <- summarizeOverlaps(ann[1], galp, mode)
    checkIdentical(1L, .getCounts(res_ga))
    checkIdentical(1L, .getCounts(res_galp))
    res_ga <- summarizeOverlaps(ann[2], ga, mode)
    res_galp <- summarizeOverlaps(ann[2], galp, mode)
    checkIdentical(1L, .getCounts(res_ga))
    checkIdentical(1L, .getCounts(res_galp))
    res_ga <- summarizeOverlaps(ann[3], ga, mode)
    res_galp <- summarizeOverlaps(ann[3], galp, mode)
    checkIdentical(1L, .getCounts(res_ga))
    checkIdentical(1L, .getCounts(res_galp))
    ## > 1 feature
    res_ga <- summarizeOverlaps(ann, ga, mode)
    res_galp <- summarizeOverlaps(ann, galp, mode)
    checkIdentical(c(0L, 0L, 0L), .getCounts(res_ga))
    checkIdentical(c(0L, 0L, 0L), .getCounts(res_galp))
    ann <- GRanges("chr1", IRanges(c(1, 2), c(23, 23)), "+")
    res_ga <- summarizeOverlaps(ann, ga, mode)
    res_galp <- summarizeOverlaps(ann, galp, mode)
    checkIdentical(c(0L, 0L), .getCounts(res_ga))
    checkIdentical(c(0L, 0L), .getCounts(res_galp))
    ann <- GRanges("chr1", IRanges(c(1, 21), c(23, 30)), "+")
    res_ga <- summarizeOverlaps(ann, ga, mode)
    res_galp <- summarizeOverlaps(ann, galp, mode)
    checkIdentical(c(0L, 0L), .getCounts(res_ga))
    checkIdentical(c(0L, 0L), .getCounts(res_galp))
    ann <- GRanges("chr1", IRanges(c(1, 1), c(23, 21)), "+")
    res_ga <- summarizeOverlaps(ann, ga, mode)
    res_galp <- summarizeOverlaps(ann, galp, mode)
    checkIdentical(c(1L, 0L), .getCounts(res_ga))
    checkIdentical(c(1L, 0L), .getCounts(res_galp))
}

test_summarizeOverlaps_features <- function()
{
    fts <- GRanges(c(rep("chr1", 7), rep("chr2", 4)), 
        IRanges(c(1000, 3000, 3600, 4000, 4000, 5000, 
            5400, 2000, 3000, 7000, 7500), 
            width = c(500, 500, 300, 500, 900, 500, 500, 
            900, 500, 600, 300)), "+",
        group = c("A", "B", "C", "C", "D", "D", "E", "F", "G", "H", "H"))
    ftslst <- split(fts, mcols(fts)[["group"]])

    rds <- GappedAlignments(c(rep(c("chr1", "chr2"), 3), "chr1"),
        as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
        c("500M", "100M", "300M", "500M", "300M", "50M200N50M", "50M150N50M"),
        strand(rep("+", 7)))

    ## GRanges
    res <- summarizeOverlaps(fts, rds, "Union")
    checkIdentical(c(1L, rep(0L, 6), 1L, 1L, 0L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(fts, rds, mode="IntersectionStrict")
    checkIdentical(c(rep(0L, 5), 1L, 0L, rep(1L, 3), 0L), .getCounts(res)) 
    res <- summarizeOverlaps(fts, rds, mode="IntersectionNotEmpty")
    checkIdentical(c(1L, rep(0L, 4), 1L, 0L, rep(1L, 3), 0L), .getCounts(res)) 

    ## GRangesList
    res <- summarizeOverlaps(ftslst, rds, "Union")
    checkIdentical(c(1L, rep(0L, 4), rep(1L, 3)), .getCounts(res)) 
    res <- summarizeOverlaps(ftslst, rds, "IntersectionStrict")
    checkIdentical(c(0L, 0L, 0L, 1L, 0L, rep(1L, 3)), .getCounts(res)) 
    res <- summarizeOverlaps(ftslst, rds, "IntersectionNotEmpty")
    checkIdentical(c(1L, 0L, 0L, 1L, 0L, rep(1L, 3)), .getCounts(res)) 
}
