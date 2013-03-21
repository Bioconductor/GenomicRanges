quiet <- suppressWarnings

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


test_GAlignmentsList_construction <- function() {
    checkTrue(validObject(GAlignmentsList()))
    checkTrue(validObject(new("GAlignmentsList")))
    checkTrue(validObject(GAlignmentsList(noGaps, Gaps)))
    checkTrue(validObject(GAlignmentsList(GappedAlignments())))
    checkTrue(validObject(GAlignmentsList(a=GappedAlignments())))
    checkException(GAlignmentsList(GRanges()), silent = TRUE)
}

test_GAlignmentsList_coercion <- function() {
    GAList <- GAlignmentsList(a=noGaps[seqnames(noGaps) == "chr3"], 
                              b=Gaps[seqnames(Gaps) == "chr4"])
    ## Lists
    rgl <- rglist(GAList)
    grl <- grglist(GAList)
    checkIdentical(length(rgl), length(GAList))
    checkIdentical(length(grl), length(GAList))
    checkIdentical(elementLengths(rgl), elementLengths(grl))
    checkIdentical(elementLengths(rgl)[1], elementLengths(GAList)[1])
    checkIdentical(elementLengths(grl)[1], elementLengths(GAList)[1])

    ## Ranges
    checkIdentical(length(ranges(GAList)), 
                   length(ranges(GAList[1])) +
                   length(ranges(GAList[2])))
    checkIdentical(length(quiet(granges(GAList))), 
                   length(quiet(granges(GAList[1]))) +
                   length(quiet(granges(GAList[2]))))
    checkIdentical(length(granges(GAList, ignore.strand=TRUE)), 
                   length(granges(GAList[1], ignore.strand=TRUE)) +
                   length(granges(GAList[2], ignore.strand=TRUE)))

    gr <- granges(GAList, ignore.strand=TRUE)
    ir <- ranges(GAList)
    checkIdentical(length(gr), length(ir))
    gr <- quiet(granges(GAList, ignore.strand=FALSE))
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
    GAList <- GAlignmentsList(a=noGaps[1:2], b=Gaps[1:2])
    checkTrue(all.equal(as.data.frame(GAList), df))
}

test_GAlignmentsList_accessors <- function() {
    GAList <- GAlignmentsList(noGaps, Gaps) 
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
    GAList <- GAlignmentsList(noGaps, Gaps)
    score <- 1:length(togroup(GAList))
    meta <- DataFrame(score=score, more=score+10) 
    mcols(GAList@unlistData) <- meta

    ## 'c' 
    checkIdentical(GAlignmentsList(), 
                   c(GAlignmentsList(), GAlignmentsList()))
    checkIdentical(GAlignmentsList(noGaps, Gaps), 
                   quiet(c(GAlignmentsList(noGaps), GAlignmentsList(Gaps))))

    ## '['
    checkIdentical(GAList, GAList[])
    checkIdentical(GAList, GAList[Rle(TRUE)])
    checkIdentical(GAList[c(TRUE, FALSE),], GAList[1])
}
