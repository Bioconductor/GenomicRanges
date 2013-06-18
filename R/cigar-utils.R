### =========================================================================
### Low-level CIGAR utilities
### -------------------------------------------------------------------------


### See p. 4 of the SAM Spec v1.4 at http://samtools.sourceforge.net/ for the
### list of CIGAR operations and their meanings.
CIGAR_OPS <- c("M", "I", "D", "N", "S", "H", "P", "=", "X")

.normarg_cigar <- function(cigar)
{
    if (is.factor(cigar))
        cigar <- as.character(cigar)
    if (!is.character(cigar))
        stop("'cigar' must be a character vector or factor")
    cigar
}

.normarg_flag <- function(flag, cigar)
{
    if (!is.null(flag)) {
        if (!is.numeric(flag))
            stop("'flag' must be NULL or a vector of integers")
        if (!is.integer(flag))
            flag <- as.integer(flag)
        if (length(cigar) != length(flag))
            stop("'cigar' and 'flag' must have the same length")
    }
    flag
}

.normarg_ops <- function(ops)
{
    if (is.null(ops))
        return(ops)
    if (!is.character(ops))
        stop("'ops' must be a character vector")
    if (any(is.na(ops)))
        stop("'ops' cannot contain NAs")
    if (length(ops) == 1L) {
        ops <- strsplit(ops, NULL, fixed=TRUE)[[1L]]
    } else if (any(nchar(ops) != 1L)) {
        stop("when 'length(ops) != 1', all its elements ",
             "must be single letters")
    }
    if (anyDuplicated(ops))
        stop("'ops' cannot contain duplicated letters")
    if (!all(ops %in% CIGAR_OPS))
        stop("'ops' contains invalid CIGAR operations")
    ops
}
 
.normarg_pos <- function(pos, cigar)
{
    if (!is.numeric(pos))
        stop("'pos' must be a vector of integers")
    if (!is.integer(pos))
        pos <- as.integer(pos)
    if (length(pos) != 1L && length(pos) != length(cigar))
        stop("'pos' must have length 1 or the same length as 'cigar'")
    pos
}

validCigar <- function(cigar)
{
    cigar <- .normarg_cigar(cigar)
    .Call2("valid_cigar", cigar, 0L, PACKAGE="GenomicRanges")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Transform CIGARs into other useful representations
###

explodeCigarOps <- function(cigar)
{
    cigar <- .normarg_cigar(cigar)
    .Call2("explode_cigar_ops", cigar, PACKAGE="GenomicRanges")
}

explodeCigarOpLengths <- function(cigar)
{
    cigar <- .normarg_cigar(cigar)
    .Call2("explode_cigar_op_lengths", cigar, PACKAGE="GenomicRanges")
}

cigarToRleList <- function(cigar)
{
    cigar_ops <- explodeCigarOps(cigar)
    cigar_op_lengths <- explodeCigarOpLengths(cigar)
    if (length(cigar) == 0L) {
        unlisted_cigar_ops <- character(0)
        unlisted_cigar_op_lengths <- integer(0)
    } else {
        unlisted_cigar_ops <- unlist(cigar_ops, use.names=FALSE)
        unlisted_cigar_op_lengths <- unlist(cigar_op_lengths, use.names=FALSE)
    }

    ## Prepare 'ans_flesh'.
    ans_flesh <- Rle(unlisted_cigar_ops, unlisted_cigar_op_lengths)

    ## Prepare 'ans_skeleton'.
    nops_per_cigar <- elementLengths(cigar_op_lengths)
    ans_breakpoints <- cumsum(unlisted_cigar_op_lengths)[cumsum(nops_per_cigar)]
    ans_skeleton <- PartitioningByEnd(ans_breakpoints)

    ## Relist.
    relist(ans_flesh, ans_skeleton)
}

splitCigar <- function(cigar)
{
    cigar <- .normarg_cigar(cigar)
    .Call2("split_cigar", cigar, PACKAGE="GenomicRanges")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### From CIGARs to ranges
###

cigarRangesOnReferenceSpace <- function(cigar, flag=NULL, ops=CIGAR_OPS,
                                        pos=1L, f=NULL,
                                        drop.empty.ranges=FALSE,
                                        reduce.ranges=FALSE,
                                        with.ops=FALSE)
{
    cigar <- .normarg_cigar(cigar)
    flag <- .normarg_flag(flag, cigar)
    ops <- .normarg_ops(ops)
    pos <- .normarg_pos(pos, cigar)
    if (!is.null(f)) {
        if (!is.factor(f))
            stop("'f' must be NULL or a factor")
        if (length(f) != length(cigar))
            stop("'f' must have the same length as 'cigar'")
    }
    if (!isTRUEorFALSE(drop.empty.ranges))
        stop("'drop.empty.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(reduce.ranges))
        stop("'reduce.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.ops))
        stop("'with.ops' must be TRUE or FALSE")
    C_ans <- .Call2("cigar_ranges",
                    cigar, flag, ops, 4L, pos, f,
                    drop.empty.ranges, reduce.ranges, with.ops,
                    PACKAGE="GenomicRanges")
    if (is.null(f))
        return(C_ans)
    compress <- length(C_ans) >= 200L
    IRangesList(C_ans, compress=compress)
}

cigarRangesOnQuerySpace <- function(cigar, flag=NULL, ops=CIGAR_OPS,
                                    drop.empty.ranges=FALSE,
                                    reduce.ranges=FALSE,
                                    with.ops=FALSE)
{
    cigar <- .normarg_cigar(cigar)
    flag <- .normarg_flag(flag, cigar)
    ops <- .normarg_ops(ops)
    if (!isTRUEorFALSE(drop.empty.ranges))
        stop("'drop.empty.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(reduce.ranges))
        stop("'reduce.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.ops))
        stop("'with.ops' must be TRUE or FALSE")
    .Call2("cigar_ranges",
           cigar, flag, ops, 1L, 1L, NULL,
           drop.empty.ranges, reduce.ranges, with.ops,
           PACKAGE="GenomicRanges")
}

cigarRangesOnPairwiseSpace <- function(cigar, flag=NULL, ops=CIGAR_OPS,
                                       drop.empty.ranges=FALSE,
                                       reduce.ranges=FALSE,
                                       with.ops=FALSE)
{
    cigar <- .normarg_cigar(cigar)
    flag <- .normarg_flag(flag, cigar)
    ops <- .normarg_ops(ops)
    if (!isTRUEorFALSE(drop.empty.ranges))
        stop("'drop.empty.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(reduce.ranges))
        stop("'reduce.ranges' must be TRUE or FALSE")
    if (!isTRUEorFALSE(with.ops))
        stop("'with.ops' must be TRUE or FALSE")
    .Call2("cigar_ranges",
           cigar, flag, ops, 3L, 1L, NULL,
           drop.empty.ranges, reduce.ranges, with.ops,
           PACKAGE="GenomicRanges")
}

### 2 wrappers to cigarRangesOnReferenceSpace() with ugly names. Used
### internally in the GenomicRanges package for turning a GAlignments object
### into a GRangesList object and/or for extracting the ranges that generate
### coverage on the reference.
cigarToIRangesListByAlignment <- function(cigar, pos=1L, flag=NULL,
                                          drop.D.ranges=FALSE,
                                          drop.empty.ranges=FALSE,
                                          reduce.ranges=TRUE)
{
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    ## It doesn't really make sense to include "I" operations here since they
    ## don't generate coverage on the reference (they always produce zero-width
    ## ranges on the reference). Anyway, this is how
    ## cigarToIRangesListByAlignment() has been behaving since the beginning.
    if (drop.D.ranges) {
        ops <- c("M", "=", "X", "I")
    } else {
        ops <- c("M", "=", "X", "I", "D")
    }
    cigarRangesOnReferenceSpace(cigar, flag=flag, ops=ops, pos=pos,
                                drop.empty.ranges=drop.empty.ranges,
                                reduce.ranges=reduce.ranges)
}

cigarToIRangesListByRName <- function(cigar, rname, pos=1L, flag=NULL,
                                      drop.D.ranges=FALSE,
                                      drop.empty.ranges=FALSE,
                                      reduce.ranges=TRUE)
{
    if (!is.factor(rname))
        stop("'rname' must be a factor")
    if (length(rname) != length(cigar))
        stop("'rname' must have the same length as 'cigar'")
    if (!isTRUEorFALSE(drop.D.ranges))
        stop("'drop.D.ranges' must be TRUE or FALSE")
    ## It doesn't really make sense to include "I" operations here since they
    ## don't generate coverage on the reference (they always produce zero-width
    ## ranges on the reference). Anyway, this is how
    ## cigarToIRangesListByRName() has been behaving since the beginning
    ## and the unit tests expect this.
    if (drop.D.ranges) {
        ops <- c("M", "=", "X", "I")
    } else {
        ops <- c("M", "=", "X", "I", "D")
    }
    cigarRangesOnReferenceSpace(cigar, flag=flag, ops=ops,
                                pos=pos, f=rname,
                                drop.empty.ranges=drop.empty.ranges,
                                reduce.ranges=reduce.ranges)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### From CIGARs to sequence lengths
###

cigarWidthOnReferenceSpace <- function(cigar)
{
    cigar <- .normarg_cigar(cigar)
    .Call2("cigar_width", cigar, 4L, PACKAGE="GenomicRanges")
}

cigarWidthOnQuerySpace <- function(cigar, before.hard.clipping=FALSE)
{
    cigar <- .normarg_cigar(cigar)
    if (!isTRUEorFALSE(before.hard.clipping))
        stop("'before.hard.clipping' must be TRUE or FALSE")
    if (before.hard.clipping) {
        space <- 0L  # QUERY_BEFORE_HARD_CLIPPING space
    } else {
        space <- 1L  # QUERY space
    }
    .Call2("cigar_width", cigar, space, PACKAGE="GenomicRanges")
}

cigarWidthOnPairwiseSpace <- function(cigar)
{
    cigar <- .normarg_cigar(cigar)
    .Call2("cigar_width", cigar, 3L, PACKAGE="GenomicRanges")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Narrow CIGARs
###

cigarNarrow <- function(cigar, start=NA, end=NA, width=NA)
{
    cigar_width <- cigarWidthOnReferenceSpace(cigar)
    cigar_ranges <- IRanges(start=rep.int(1L, length(cigar_width)),
                            width=cigar_width)
    threeranges <- threebands(cigar_ranges, start=start, end=end, width=width)
    C_ans <- .Call2("cigar_narrow",
                   cigar, width(threeranges$left), width(threeranges$right),
                   PACKAGE="GenomicRanges")
    ans <- C_ans[[1L]]
    attr(ans, "rshift") <- C_ans[[2L]]
    ans
}

cigarQNarrow <- function(cigar, start=NA, end=NA, width=NA)
{
    cigar_qwidth <- cigarWidthOnQuerySpace(cigar)
    cigar_qranges <- IRanges(start=rep.int(1L, length(cigar_qwidth)),
                             width=cigar_qwidth)
    threeranges <- threebands(cigar_qranges, start=start, end=end, width=width)
    C_ans <- .Call2("cigar_qnarrow",
                   cigar, width(threeranges$left), width(threeranges$right),
                   PACKAGE="GenomicRanges")
    ans <- C_ans[[1L]]
    attr(ans, "rshift") <- C_ans[[2L]]
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Summarize CIGARs
###

cigarOpTable <- function(cigar)
{
    cigar <- .normarg_cigar(cigar)
    .Call2("cigar_op_table", cigar, PACKAGE="GenomicRanges")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Translate coordinates between query-based and reference-based
###

queryLoc2refLoc <- function(qloc, cigar, pos=1L)
{
    stop("NOT IMPLEMENTED YET, SORRY!")
}

queryLocs2refLocs <- function(qlocs, cigar, pos=1L, flag=NULL)
{
    stop("NOT IMPLEMENTED YET, SORRY!")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Old stuff (deprecated & defunct)
###

cigarToWidth <- function(...)
{
    .Deprecated("cigarWidthOnReferenceSpace")
    cigarWidthOnReferenceSpace(...)
}

cigarToQWidth <- function(...)
{
    .Deprecated("cigarWidthOnQuerySpace")
    cigarWidthOnQuerySpace(...)
}

cigarToIRanges <- function(cigar, drop.D.ranges=FALSE,
                           drop.empty.ranges=FALSE, reduce.ranges=TRUE)
{
    .Deprecated("cigarToIRangesListByAlignment")
    cigarToIRangesListByAlignment(cigar, drop.D.ranges=drop.D.ranges,
                                  drop.empty.ranges=drop.empty.ranges,
                                  reduce.ranges=reduce.ranges)[[1L]]
}

cigarToCigarTable <- function(cigar)
{
    msg <- c("  cigarToCigarTable() is deprecated.",
             "Please use 'table(cigar)' instead.")
    .Deprecated(msg=paste0(msg, collapse=" "))
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

### summarizeCigarTable() is broken:
###
###   > summarizeCigarTable(cigarToCigarTable(c("55M", "5M3I40M")))
###   Error in data.frame(subject = subjectHits(findOverlaps(IRanges(unlist(x[["cigar"]] ==  : 
###     arguments imply differing number of rows: 0, 1
###
### In addition, what it's supposed to do and to return is not really
### documented, its implementation is kind of cryptic, and it doesn't seem
### very useful anyway (hard to know exactly without more details about what
### it does though). I doubt a lot of people have tried to use it ==> Probably
### not worth fixing/maintaining.
summarizeCigarTable <- function(x)
{
    msg <- "  summarizeCigarTable() is deprecated."
    .Deprecated(msg=msg)
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

