.getCounts <- function(res)
{
    as.integer(assays(res)$counts)
}
quiet <- suppressMessages

gr <- GRanges(c(rep("chr1", 7), rep("chr2", 4)), 
    IRanges(c(1000, 3000, 3600, 4000, 4000, 5000, 
        5400, 2000, 3000, 7000, 7500), 
        width = c(500, 500, 300, 500, 900, 500, 500, 
        900, 500, 600, 300)), "+",
    group = c("A", "B", "C", "C", "D", "D", "E", "F", "G", "H", "H"))

rds <- GAlignments(c(rep(c("chr1", "chr2"), 3), "chr1"),
    as.integer(c(1400, 2700, 3400, 7100, 4000, 3100, 5200)),
    c("500M", "100M", "300M", "500M", "300M", "50M200N50M", "50M150N50M"),
    strand(rep("+", 7)))

test_summarizeOverlaps_Union_single <- function()
{
    ## single-end no gaps
    mode <- "Union"
    ga <- GAlignments("chr1", 20L, "11M", strand("+")) 
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
    ga <- GAlignments("chr1", 1L, "10M4N11M", strand("+"))
    ga1 <- GAlignments("chr1", 1L, "10M", strand("+"))
    ga2 <- GAlignments("chr1", 15L, "11M", strand("-"))
    galp <- GAlignmentPairs(ga1, ga2, TRUE)
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
    ga <- GAlignments("chr1", 7L, "6M", strand("+")) 
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
    ga <- GAlignments("chr1", 10L, "6M4N6M", strand("+"))
    ga1 <- GAlignments("chr1", 10L, "6M", strand("+"))
    ga2 <- GAlignments("chr1", 20L, "6M", strand("-"))
    galp <- GAlignmentPairs(ga1, ga2, TRUE)
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
    ga <- GAlignments("chr1", 10L, "11M", strand("+")) 
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
    ga <- GAlignments("chr1", 23L, "5M", strand("*"))
    res <- summarizeOverlaps(ann, ga, mode)
    checkIdentical(c(0L, 0L, 0L), .getCounts(res))
    strand(ga) <- "-"
    res <- summarizeOverlaps(ann, ga, mode)
    checkIdentical(c(0L, 0L, 1L), .getCounts(res))
    ga <- GAlignments("chr1", 28L, "5M", strand("+"))
    res <- summarizeOverlaps(ann, ga, mode)
    checkIdentical(c(1L, 0L, 0L), .getCounts(res))
}

test_summarizeOverlaps_IntersectionNotEmpty_paired <- function()
{
    ## single-end with gaps behave like paired-end
    mode <- "IntersectionNotEmpty"
    ga <- GAlignments("chr1", 10L, "6M4N6M", strand("+"))
    ga1 <- GAlignments("chr1", 10L, "6M", strand("+"))
    ga2 <- GAlignments("chr1", 20L, "6M", strand("-"))
    galp <- GAlignmentPairs(ga1, ga2, TRUE)
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

test_summarizeOverlaps_inter.feature_GRanges <- function()
{
    ## rows 5,6,7 from figure in vignette
    ft <- gr[10:11]
    rd <- GAlignments(rep("chr2", 3), as.integer(c(7100, 7100, 7500)), 
                      c("300M", "500M", "50M"), strand(rep("+", 3)))
    mode <- "Union"
    res <- summarizeOverlaps(ft, rd[1], mode, inter.feature=TRUE)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[1], mode, inter.feature=FALSE)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[2], mode, inter.feature=TRUE)
    checkIdentical(c(0L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[2], mode, inter.feature=FALSE)
    checkIdentical(c(1L, 1L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[3], mode, inter.feature=TRUE)
    checkIdentical(c(0L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[3], mode, inter.feature=FALSE)
    checkIdentical(c(1L, 1L), .getCounts(res)) 

    mode <- "IntersectionStrict"
    res <- summarizeOverlaps(ft, rd[1], mode, inter.feature=TRUE)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[1], mode, inter.feature=FALSE)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[2], mode, inter.feature=TRUE)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[2], mode, inter.feature=FALSE)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[3], mode, inter.feature=TRUE)
    checkIdentical(c(0L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[3], mode, inter.feature=FALSE)
    checkIdentical(c(1L, 1L), .getCounts(res)) 

    mode <- "IntersectionNotEmpty"
    res <- summarizeOverlaps(ft, rd[1], mode, inter.feature=TRUE)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[1], mode, inter.feature=FALSE)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[2], mode, inter.feature=TRUE)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[2], mode, inter.feature=FALSE)
    checkIdentical(c(1L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[3], mode, inter.feature=TRUE)
    checkIdentical(c(0L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd[3], mode, inter.feature=FALSE)
    checkIdentical(c(0L, 0L), .getCounts(res))
    ## read spans both features 
    rd <- GAlignments("chr2", 7000L, "750M", strand("+"))
    res <- summarizeOverlaps(ft, rd, mode, inter.feature=TRUE)
    checkIdentical(c(0L, 0L), .getCounts(res)) 
    res <- summarizeOverlaps(ft, rd, mode, inter.feature=FALSE)
    checkIdentical(c(1L, 1L), .getCounts(res)) 
}

test_summarizeOverlaps_inter.feature_GRangesList <- function()
{
    grl <- split(gr, mcols(gr)[["group"]])
    mode <- "Union"
    res <- quiet(summarizeOverlaps(grl, rds, mode))
    checkIdentical(c(1L, rep(0L, 4), rep(1L, 3)), .getCounts(res)) 
    res <- quiet(summarizeOverlaps(grl, rds, mode, inter.feature=FALSE))
    checkIdentical(c(1L, 1L, 2L, 2L, rep(1L, 4)), .getCounts(res)) 
    co <- countOverlaps(grl, rds, type="any")
    checkIdentical(unname(co), .getCounts(res))

    mode <- "IntersectionStrict"
    res <- quiet(summarizeOverlaps(grl, rds, mode))
    checkIdentical(c(0L, 0L, 0L, 1L, 0L, rep(1L, 3)), .getCounts(res)) 
    res <- quiet(summarizeOverlaps(grl, rds, mode, inter.feature=FALSE))
    checkIdentical(c(0L, 0L, 1L, 2L, 0L, rep(1L, 3)), .getCounts(res))
    co <- countSubjectHits(findOverlaps(rds, grl, type="within"))
    checkIdentical(unname(co), .getCounts(res))
 
    mode <- "IntersectionNotEmpty"
    rd <-  rds[c(5,7)]
    ft <- grl[3:5]
    res <- quiet(summarizeOverlaps(ft, rd, mode))
    checkIdentical(c(0L, 1L, 0L), .getCounts(res)) 
    res <- quiet(summarizeOverlaps(ft, rd, mode, inter.feature=FALSE))
    checkIdentical(c(0L, 1L, 0L), .getCounts(res)) 
}

