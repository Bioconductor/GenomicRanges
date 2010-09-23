### =========================================================================
### The GenomicRanges interface
### -------------------------------------------------------------------------
###

setClass("GenomicRanges", representation("VIRTUAL"))

### The code in this file will work out-of-the-box on 'x' as long as
### seqnames(x), ranges(x), strand(x), seqlengths(x) and isCircular()
### are defined.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GenomicRanges.length <- function(x)
{
    n <- length(seqnames(x))
    if ((length(ranges(x)) != n) || (length(strand(x)) != n) ||
        (nrow(elementMetadata(x)) != n))
        "slot lengths are not all equal"
    else
        NULL
}

.valid.GenomicRanges.seqnames <- function(x)
{
    if (!is.factor(runValue(seqnames(x))))
        "slot 'seqnames' should be a 'factor' Rle"
    else if (IRanges:::anyMissing(runValue(seqnames(x))))
        "slot 'seqnames' contains missing values"
    else
        NULL
}

.valid.GenomicRanges.strand <- function(x)
{
    if (!is.factor(runValue(strand(x))) ||
        !identical(levels(runValue(strand(x))), levels(strand())))
        paste("slot 'strand' should be a 'factor' Rle with levels c(",
              paste('"', levels(strand()), '"', sep = "", collapse = ", "),
                    ")", sep = "")
    else if (IRanges:::anyMissing(runValue(strand(x))))
        "slot 'strand' contains missing values"
    else
        NULL
}

.valid.GenomicRanges.seqlengths <- function(x)
{
    msg <- NULL
    if (length(x) > 0 && is.null(names(seqlengths(x))))
        msg <- "slot 'seqlengths' is unnamed"
    if (!setequal(names(seqlengths(x)), levels(seqnames(x))))
        msg <-
          c(msg, "slot 'seqlengths' names do not match 'levels(seqnames)'")
    if (IRanges:::anyMissing(seqlengths(x))) {
        if (!all(is.na(seqlengths(x))))
            msg <-
              c(msg,
                "slot 'seqlengths' cannot mix NAs and non-NA integers")
    } else {
        if (any(seqlengths(x) < 0L, na.rm = TRUE))
            msg <- c(msg, "slot 'seqlengths' contains negative values")
        ## TODO: Loosen the check below for circular sequences.
        seqnames <- seqnames(x)
        runValue(seqnames) <- runValue(seqnames)[drop=TRUE]
        minStarts <- IRanges:::.tapplyDefault(start(x), seqnames, min)
        maxEnds <- IRanges:::.tapplyDefault(end(x), seqnames, max)
        if (any(minStarts < 1L) || any(maxEnds > seqlengths(x)[names(maxEnds)]))
            msg <-
              c(msg, "slot 'ranges' contains values outside of sequence bounds")
    }
    msg
}

.valid.GenomicRanges.isCircular <- function(x)
{
    msg <- NULL
    if (length(x) > 0 && is.null(names(isCircular(x))))
        msg <- "slot 'isCircular' is unnamed"
    if (!setequal(names(isCircular(x)), levels(seqnames(x))))
        msg <-
          c(msg, "slot 'isCircular' names do not match 'levels(seqnames)'")
    if (IRanges:::anyMissing(isCircular(x))) {
        if (!all(is.na(isCircular(x))))
            msg <-
              c(msg,
                "slot 'isCircular' cannot mix NAs and non-NA integers")
    }
    msg
}

.valid.GenomicRanges.elementMetadata <- function(x)
{
    msg <- NULL
    ## NOTE: This list is also included in the man page for GRanges objects.
    ## Keep the 2 lists in sync!
    INVALID.COLNAMES <- c("seqnames", "ranges", "strand",
                          "seqlengths", "isCircular",
                          "start", "end", "width", "element")
    if (any(INVALID.COLNAMES %in% colnames(elementMetadata(x))))
        msg <-
          paste("slot 'elementMetadata' cannot use",
                paste("\"", INVALID.COLNAMES, "\"", sep="", collapse=", "),
                "as column names")
    msg
}

.valid.GenomicRanges <- function(x)
{
    c(.valid.GenomicRanges.length(x),
      .valid.GenomicRanges.seqnames(x),
      .valid.GenomicRanges.strand(x),
      .valid.GenomicRanges.seqlengths(x),
      .valid.GenomicRanges.isCircular(x),
      .valid.GenomicRanges.elementMetadata(x))
}

setValidity2("GenomicRanges", .valid.GenomicRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("GenomicRanges", "RangedData",
      function(from)
      {
        rd <- RangedData(ranges(from), strand = strand(from),
                         elementMetadata(from), space = seqnames(from))
        elementMetadata(ranges(rd)) <- DataFrame(
                                         seqlengths = seqlengths(from),
                                         isCircular = isCircular(from))
        rd
      }
)

setAs("GenomicRanges", "RangesList",
      function(from)
      {
        rl <- split(ranges(from), seqnames(from))
        emd <- split(DataFrame(strand = strand(from), elementMetadata(from)),
                     seqnames(from))
        rl <- mseqapply(function(ranges, metadata) {
          elementMetadata(ranges) <- metadata
          ranges
        }, rl, emd)
        elementMetadata(rl) <- DataFrame(seqlengths = seqlengths(from),
                                         isCircular = isCircular(from))
        rl
      }
      )

setMethod("as.data.frame", "GenomicRanges",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        ranges <- ranges(x)
        if (missing(row.names))
            row.names <- names(x)
        if (!is.null(names(x)))
            names(x) <- NULL
        data.frame(seqnames = as.factor(seqnames(x)),
                   start = start(x),
                   end = end(x),
                   width = width(x),
                   strand = as.factor(strand(x)),
                   as.data.frame(elementMetadata(x)),
                   row.names = row.names,
                   stringsAsFactors = FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("elementMetadata", "GenomicRanges",
    function(x, ...)
    {
        ans <- callNextMethod()
        if (!is.null(names(x)))
            rownames(ans) <- names(x)
        ans
    }
)

setMethod("names", "GenomicRanges", function(x) names(ranges(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ranges methods.
###

setMethod("start", "GenomicRanges", function(x, ...) start(ranges(x)))
setMethod("end", "GenomicRanges", function(x, ...) end(ranges(x)))
setMethod("width", "GenomicRanges", function(x) width(ranges(x)))

setMethod("coverage", "GenomicRanges",
    function(x, shift = list(0L), width = as.list(seqlengths(x)),
             weight = list(1L))
    {
        fixArg <- function(arg, argname, uniqueSeqnames) {
            k <- length(uniqueSeqnames)
            if (!is.list(arg) || is(arg, "IntegerList"))
                stop("'", argname, "' must be a list")
            makeNULL <-
              unlist(lapply(arg, function(y) length(y) == 1 && is.na(y)),
                     use.names=FALSE)
            if (any(makeNULL))
                arg[makeNULL] <- rep(list(NULL), sum(makeNULL))
            if (length(arg) < k)
                arg <- rep(arg, length.out = k)
            if (is.null(names(arg)))
                names(arg) <- uniqueSeqnames
            if (!all(uniqueSeqnames %in% names(arg)))
                stop("some seqnames missing from names(", argname, ")")
            arg
        }
        uniqueSeqnames <- levels(seqnames(x))
        shift <- fixArg(shift, "shift", uniqueSeqnames)
        width <- fixArg(width, "width", uniqueSeqnames)
        weight <- fixArg(weight, "weight", uniqueSeqnames)
        xSplitRanges <- splitRanges(seqnames(x))
        xRanges <- unname(ranges(x))
        IRanges:::newSimpleList("SimpleRleList",
                      lapply(structure(uniqueSeqnames, names = uniqueSeqnames),
                             function(i) {
                                 coverage(seqselect(xRanges, xSplitRanges[[i]]),
                                          shift = shift[[i]], width = width[[i]],
                                          weight = weight[[i]])
                             }))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sequence methods.
###

setMethod("length", "GenomicRanges", function(x) length(seqnames(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

.showSeqlengths <-
function(object)
{
    seqlens <- seqlengths(object)
    nseq <- length(seqlens)
    halfWidth <- getOption("width") %/% 2L
    first <- max(1L, halfWidth)
    showMatrix <-
      rbind(as.character(head(names(seqlens), first)),
            as.character(head(seqlens, first)))
    if (nseq > first) {
        last <- min(nseq - first, halfWidth)
        showMatrix <-
          cbind(showMatrix,
                rbind(as.character(tail(names(seqlens), last)),
                      as.character(tail(seqlens, last))))
    }
    showMatrix <- format(showMatrix, justify = "right")
    cat("\nseqlengths\n")
    cat(IRanges:::labeledLine("", showMatrix[1L,], count = FALSE,
        labelSep = ""))
    cat(IRanges:::labeledLine("", showMatrix[2L,], count = FALSE,
        labelSep = ""))
}

showGenomicRanges <-
function(object, print.seqlengths = FALSE)
{
    lo <- length(object)
    nc <- ncol(elementMetadata(object))
    cat(class(object), " with ",
        lo, ifelse(lo == 1, " range and ", " ranges and "),
        nc, ifelse(nc == 1, " elementMetadata value\n",
                   " elementMetadata values\n"),
            sep = "")
    if (lo == 0) {
        out <-
          matrix(nrow = 0L, ncol = 4L + nc,
                 dimnames = list(NULL,
                                 c("seqnames", "ranges", "strand", "|",
                                   colnames(elementMetadata(object)))))
    } else {
        nms <- names(object)
        if (lo < 20) {
            out <-
              cbind(seqnames = as.character(seqnames(object)),
                    ranges = IRanges:::showAsCell(ranges(object)),
                    strand = as.character(strand(object)),
                    "|" = rep.int("|", lo))
            if (nc > 0)
                out <-
                  cbind(out,
                        as.matrix(format(do.call(data.frame,
                                                 lapply(elementMetadata(object),
                                                        IRanges:::showAsCell)))))
            if (is.null(nms))
                rownames(out) <-
                  format(paste("[", seq_len(lo), "]", sep = ""),
                         justify = "right")
            else
                rownames(out) <- format(nms, justify = "right")
            classinfo <-
              matrix(c("<Rle>", "<IRanges>", "<Rle>", "|",
                       unlist(lapply(elementMetadata(object), function(x)
                                     paste("<", class(x), ">", sep = "")),
                              use.names = FALSE)), nrow = 1,
                     dimnames = list("", colnames(out)))
        } else {
            top <- object[1:9]
            bottom <- object[(lo-8L):lo]
            out <-
              rbind(cbind(seqnames = as.character(seqnames(top)),
                          ranges = IRanges:::showAsCell(ranges(top)),
                          strand = as.character(strand(top)),
                          "|" = rep.int("|", 9)),
                    rbind(rep.int("...", 4)),
                    cbind(seqnames = as.character(seqnames(bottom)),
                          ranges = IRanges:::showAsCell(ranges(bottom)),
                          strand = as.character(strand(bottom)),
                          "|" = rep.int("|", 9)))
            if (nc > 0)
                out <-
                  cbind(out,
                        rbind(as.matrix(format(do.call(data.frame,
                                                       lapply(elementMetadata(top),
                                                              IRanges:::showAsCell)))),
                                        rbind(rep.int("...", nc)),
                                        rbind(as.matrix(format(do.call(data.frame,
                                                                        lapply(elementMetadata(bottom),
                                                                               IRanges:::showAsCell)))))))
            if (is.null(nms)) {
                rownames(out) <-
                  format(c(paste("[", 1:9, "]", sep = ""), "...",
                           paste("[", (lo-8L):lo, "]", sep = "")),
                         justify = "right")
            } else {
                rownames(out) <-
                  format(c(head(nms, 9), "...", tail(nms, 9)),
                         justify = "right")
            }
            classinfo <-
              matrix(c("<Rle>", "<IRanges>", "<Rle>", "|",
                       unlist(lapply(elementMetadata(top), function(x)
                                     paste("<", class(x), ">", sep = "")),
                              use.names = FALSE)), nrow = 1,
                     dimnames = list("", colnames(out)))
        }
        out <- rbind(classinfo, out)
    }
    print(out, quote = FALSE, right = TRUE)
    if (print.seqlengths) {
        .showSeqlengths(object)
    }
}

setMethod("show", "GenomicRanges",
          function(object) showGenomicRanges(object, print.seqlengths = TRUE))

