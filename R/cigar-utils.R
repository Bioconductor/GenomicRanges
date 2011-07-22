###

validCigar <- function(cigar)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    .Call2("valid_cigar", cigar, 0L, PACKAGE="GenomicRanges")
}

cigarOpTable <- function(cigar)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    .Call2("cigar_op_table", cigar, PACKAGE="GenomicRanges")
}

cigarToQWidth <- function(cigar, before.hard.clipping=FALSE)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    if (!isTRUEorFALSE(before.hard.clipping))
        stop("'before.hard.clipping' must be TRUE or FALSE")
    .Call2("cigar_to_qwidth",
          cigar, before.hard.clipping,
          PACKAGE="GenomicRanges")
}

cigarToWidth <- function(cigar)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    .Call2("cigar_to_width", cigar, PACKAGE="GenomicRanges")
}

cigarQNarrow <- function(cigar, start=NA, end=NA, width=NA)
{
    threeranges <- threebands(successiveIRanges(cigarToQWidth(cigar)),
                              start=start, end=end, width=width)
    C_ans <- .Call2("cigar_qnarrow",
                   cigar, width(threeranges$left), width(threeranges$right),
                   PACKAGE="GenomicRanges")
    ans <- C_ans[[1L]]
    attr(ans, "rshift") <- C_ans[[2L]]
    ans
}

cigarNarrow <- function(cigar, start=NA, end=NA, width=NA)
{
    threeranges <- threebands(successiveIRanges(cigarToWidth(cigar)),
                              start=start, end=end, width=width)
    C_ans <- .Call2("cigar_narrow",
                   cigar, width(threeranges$left), width(threeranges$right),
                   PACKAGE="GenomicRanges")
    ans <- C_ans[[1L]]
    attr(ans, "rshift") <- C_ans[[2L]]
    ans
}

cigarToIRanges <- function(cigar, drop.D.ranges=FALSE, merge.ranges=TRUE)
{
    if (is.factor(cigar) && is.character(levels(cigar)))
        cigar <- as.vector(cigar)
    if (!isSingleString(cigar))
        stop("'cigar' must be a single string")
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(merge.ranges))
        stop("'merge.ranges' must be TRUE or FALSE")
    .Call2("cigar_to_IRanges", cigar, drop.D.ranges, merge.ranges,
          PACKAGE="GenomicRanges")
}

cigarToIRangesListByAlignment <-
function(cigar, pos, flag=NULL, drop.D.ranges=FALSE)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    if (!is.numeric(pos))
        stop("'pos' must be a vector of integers")
    if (!is.integer(pos))
        pos <- as.integer(pos)
    if (length(cigar) != length(pos))
        stop("'cigar' and 'pos' must have the same length")
    if (!is.null(flag)) {
        if (!is.numeric(flag))
            stop("'flag' must be NULL or a vector of integers")
        if (!is.integer(flag))
            flag <- as.integer(flag)
        if (length(cigar) != length(flag))
            stop("'cigar' and 'flag' must have the same length")
    }
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    .Call2("cigar_to_list_of_IRanges_by_alignment",
          cigar, pos, flag, drop.D.ranges, PACKAGE="GenomicRanges")
}

cigarToIRangesListByRName <-
function(cigar, rname, pos, flag=NULL, drop.D.ranges=FALSE, merge.ranges=TRUE)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    if (!is.factor(rname) || !is.character(levels(rname))) {
        if (!is.character(rname))
            stop("'rname' must be a character vector/factor")
        rname <- as.factor(rname)
    }
    if (!is.numeric(pos))
        stop("'pos' must be a vector of integers")
    if (!is.integer(pos))
        pos <- as.integer(pos)
    if (length(cigar) != length(rname) || length(cigar) != length(pos))
        stop("'cigar', 'rname' and 'pos' must have the same length")
    if (!is.null(flag)) {
        if (!is.numeric(flag))
            stop("'flag' must be NULL or a vector of integers")
        if (!is.integer(flag))
            flag <- as.integer(flag)
        if (length(cigar) != length(flag))
            stop("'cigar' and 'flag' must have the same length")
    }
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(merge.ranges))
        stop("'merge.ranges' must be TRUE or FALSE")
    C_ans <- .Call2("cigar_to_list_of_IRanges_by_rname",
                   cigar, rname, pos, flag, drop.D.ranges, merge.ranges,
                   PACKAGE="GenomicRanges")
    if (length(C_ans) < 200L)
        IRangesList(C_ans, compress=FALSE)
    else
        IRangesList(C_ans, compress=TRUE)
}

queryLoc2refLoc <- function(qloc, cigar, pos=1)
{
    stop("NOT IMPLEMENTED YET, SORRY!")
}

queryLocs2refLocs <- function(qlocs, cigar, pos, flag=NULL)
{
    stop("NOT IMPLEMENTED YET, SORRY!")
}

splitCigar <- function(cigar)
{
    if (!is.character(cigar)) {
        if (!is.factor(cigar) || !is.character(levels(cigar)))
            stop("'cigar' must be a character vector/factor")
        cigar <- as.vector(cigar)
    }
    .Call2("split_cigar", cigar, PACKAGE="GenomicRanges")
}

cigarToRleList <- function(cigar)
{
    splitCigars <- splitCigar(cigar)
    splitRawValues <- unname(lapply(splitCigars, "[[", 1L))
    splitLengths <- unname(lapply(splitCigars, "[[", 2L))
    IRanges:::newCompressedList("CompressedRleList",
                           Rle(rawToChar(unlist(splitRawValues), multiple=TRUE),
                               unlist(splitLengths)),
                           cumsum(unlist(lapply(splitLengths, sum))))
}

cigarToCigarTable <- function(cigar) {
    if (is.character(cigar))
        cigar <- factor(cigar)
    else if (!is.factor(cigar))
        stop("'cigar' must be a character vector/factor")
    basicTable <- table(cigar)
    tableOrder <- order(basicTable, decreasing=TRUE)
    cigar <- factor(cigar, levels = levels(cigar)[tableOrder])
    basicTable <- basicTable[tableOrder]
    DataFrame(cigar = cigarToRleList(levels(cigar)),
              count = as.integer(basicTable))
}

summarizeCigarTable <- function(x) {
    alignedCharacters <-
      table(rep.int(elementLengths(x[["cigar"]]), x[["count"]]),
            rep.int(viewSums(Views(unlist(x[["cigar"]] == "M"),
                                   as(x[["cigar"]]@partitioning, "IRanges"))) ==
                                   elementLengths(x[["cigar"]]),
                    x[["count"]]))
    tabledAlignedCharacters <- as(rev(colSums(alignedCharacters)), "integer")
    names(tabledAlignedCharacters) <-
      unname(c("TRUE" = "AllAligned",
               "FALSE" = "SomeNonAligned")[names(tabledAlignedCharacters)])

    indelHits <-
      rbind(data.frame(subject =
                       subjectHits(findOverlaps(IRanges(unlist(x[["cigar"]] == "D")),
                                                x[["cigar"]]@partitioning)),
                       type = factor("D", levels = c("D", "I"))),
            data.frame(subject =
                       subjectHits(findOverlaps(IRanges(unlist(x[["cigar"]] == "I")),
                                                x[["cigar"]]@partitioning)),
                       type = factor("I", levels = c("D", "I"))))
    tabledIndelHits <- table(indelHits[,1], indelHits[,2])
    tabledIndelHits <-
      tabledIndelHits[rep.int(seq_len(nrow(tabledIndelHits)),
                              x[["count"]][as.integer(rownames(tabledIndelHits))]),
                      , drop = FALSE]
    tabledIndelHits <-
      as(table(tabledIndelHits[,"D"], tabledIndelHits[,"I"]), "matrix")
    rownames(tabledIndelHits) <- paste("D", rownames(tabledIndelHits), sep = "")
    colnames(tabledIndelHits) <- paste("I", colnames(tabledIndelHits), sep = "")
    if (!("D0" %in% rownames(tabledIndelHits)))
        tabledIndelHits <-
          rbind("D0" = rep.int(0L, ncol(tabledIndelHits)), tabledIndelHits)
    if (!("I0" %in% colnames(tabledIndelHits)))
        tabledIndelHits <-
          cbind("I0" = rep.int(0L, nrow(tabledIndelHits)), tabledIndelHits)
    tabledIndelHits["D0", "I0"] <- nrow(x) - sum(tabledIndelHits[-1])

    list("AlignedCharacters" = tabledAlignedCharacters,
         "Indels" = tabledIndelHits)
}
