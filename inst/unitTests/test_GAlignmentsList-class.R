noGaps <- GappedAlignments(
    Rle(factor(c("chr1", "chr2", "chr1", "chr3")), 
        c(1, 3, 2, 4)), 
    pos=1:10, cigar=paste0(10:1, "M"),
    strand=Rle(strand(c("-", "+", "*", "+", "-")), 
        c(1, 2, 2, 3, 2)),
    names=head(letters, 10), score=1:10)
Gaps <- GappedAlignments(
    Rle(factor(c("chr2", "chr4")), c(3, 4)), pos=1:7,
    cigar=c("5M", "3M2N3M2N3M", "5M", "10M", "5M1N4M", "8M2N1M", "5M"), 
    strand=Rle(strand(c("-", "+")), c(4, 3)),
    names=tail(letters, 7), score=1:7)
GAList <- GAlignmentsList(a=noGaps, b=Gaps)
quiet <- suppressWarnings

test_GAlignmentsList_construction <- function() {
    checkTrue(validObject(GAlignmentsList()))
    checkTrue(validObject(new("GAlignmentsList")))
    checkTrue(validObject(GAlignmentsList(noGaps, Gaps)))
    checkTrue(validObject(GAlignmentsList(GappedAlignments())))
    checkTrue(validObject(GAlignmentsList(a=GappedAlignments())))
    checkException(GAlignmentsList(GRanges()), silent = TRUE)
}

test_GAlignmentsList_coercion <- function() {
    galist <- GAlignmentsList(a=noGaps[seqnames(noGaps) == "chr3"], 
                              b=Gaps[seqnames(Gaps) == "chr4"])
    ## Lists
    rgl <- rglist(galist)
    grl <- grglist(galist)
    checkIdentical(length(rgl), length(galist))
    checkIdentical(length(grl), length(galist))
    checkIdentical(elementLengths(rgl), elementLengths(grl))
    checkIdentical(elementLengths(rgl)[1], elementLengths(galist)[1])
    checkIdentical(elementLengths(grl)[1], elementLengths(galist)[1])

    ## Ranges
    checkIdentical(length(ranges(galist)), 
                   length(ranges(galist[1])) +
                   length(ranges(galist[2])))
    checkIdentical(length(quiet(granges(galist))), 
                   length(quiet(granges(galist[1]))) +
                   length(quiet(granges(galist[2]))))
    checkIdentical(length(granges(galist, ignore.strand=TRUE)), 
                   length(granges(galist[1], ignore.strand=TRUE)) +
                   length(granges(galist[2], ignore.strand=TRUE)))

    gr <- granges(galist, ignore.strand=TRUE)
    ir <- ranges(galist)
    checkIdentical(length(gr), length(ir))
    gr <- quiet(granges(galist, ignore.strand=FALSE))
    checkTrue(length(gr) == 4L)

    ## data.frame
    df <- data.frame(element=rep(c("a", "b"), each=2),
                     seqnames=c("chr1", rep("chr2", 3)), 
                     strand=c("-", "+", "-", "-"),
                     cigar=c("10M", "9M", "5M", "3M2N3M2N3M"),
                     qwidth=c(10, 9 , 5, 9), start=c(1, 2, 1, 2),
                     end=c(10, 10, 5, 14), width=c(10, 9, 5, 13),
                     ngap=c(0, 0, 0, 2), score=c(1, 2, 1, 2), 
                     row.names=c("a", "b", "t", "u"),
                     stringsAsFactors=FALSE)
    galist <- GAlignmentsList(a=noGaps[1:2], b=Gaps[1:2])
    checkTrue(all.equal(as.data.frame(galist), df))

    ## introns
    galist <- GAList
    grl <- introns(galist)
    checkIdentical(names(galist), names(grl))
    checkTrue(length(galist) == length(grl))
    checkTrue(length(grl[[1]]) == 0L)
    checkTrue(length(grl[[2]]) == 4L)
}

test_GAlignmentsList_accessors <- function() {
    galist <- GAlignmentsList(noGaps, Gaps) 
    target <- RleList(lapply(GAList, seqnames), compress=TRUE)
    checkIdentical(seqnames(GAList), target) 
    target <- RleList(lapply(GAList, rname), compress=TRUE)
    checkIdentical(rname(GAList), target)
    target <- CharacterList(lapply(GAList, cigar), compress=TRUE)
    checkIdentical(cigar(GAList), target) 
    target <- RleList(lapply(GAList, strand), compress=TRUE)
    checkIdentical(strand(GAList), target) 
    target <- IntegerList(lapply(GAList, width))
    checkIdentical(width(GAList), target)
    target <- SplitDataFrameList(lapply(GAList, mcols))
    checkIdentical(mcols(GAList, level="within"), target)
}

test_GAlignments_subset_combine <- function()
{
    galist <- GAList
    score <- 1:length(togroup(galist))
    meta <- DataFrame(score=score, more=score+10) 
    mcols(galist@unlistData) <- meta

    ## 'c' 
    checkIdentical(GAlignmentsList(), 
                   c(GAlignmentsList(), GAlignmentsList()))
    checkIdentical(GAlignmentsList(noGaps, Gaps), 
                   quiet(c(GAlignmentsList(noGaps), GAlignmentsList(Gaps))))

    ## '['
    checkIdentical(galist, galist[])
    checkIdentical(galist, galist[Rle(TRUE)])
    checkIdentical(galist[c(TRUE, FALSE),], galist[1])
}

test_GAlignments_qnarrow <- function()
{
    galist <- GAlignmentsList(noGaps[1:6], Gaps)
    qn <- qnarrow(galist, end=-4)
    checkIdentical(qnarrow(galist[[1]], end=-4), qn[[1]])
    checkIdentical(qnarrow(galist[[2]], end=-4), qn[[2]])
}
