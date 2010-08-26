### =========================================================================
### GRanges objects
### -------------------------------------------------------------------------
###
### Class definition

setClass("GRanges", contains = c("Sequence", "GenomicRanges"),
         representation(seqnames = "Rle",
                        ranges = "IRanges",
                        strand = "Rle",
                        seqlengths = "integer"),
         prototype(seqnames = Rle(factor()),
                   strand = Rle(strand()),
                   elementMetadata = DataFrame()))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GRanges.elementMetadata <- function(x)
{
    msg <- NULL
    if (!is.null(rownames(x@elementMetadata)))
        msg <- c(msg, "slot 'elementMetadata' cannot contain row names")
    msg
}

setValidity2("GRanges", .valid.GRanges.elementMetadata)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GRanges <-
function(seqnames = Rle(), ranges = IRanges(),
         strand = Rle("*", length(seqnames)),
         ...,
         seqlengths =
         structure(rep(NA_integer_, length(levels(seqnames))),
                   names = levels(seqnames)))
{
    if (!is(seqnames, "Rle"))
        seqnames <- Rle(seqnames)
    if (!is.factor(runValue(seqnames)))
        runValue(seqnames) <- factor(runValue(seqnames))

    if (!is(ranges, "IRanges"))
        ranges <- as(ranges, "IRanges")

    if (!is(strand, "Rle"))
        strand <- Rle(strand)
    if (!is.factor(runValue(strand)) ||
        !identical(levels(runValue(strand)), levels(strand())))
        runValue(strand) <- strand(runValue(strand))
    if (IRanges:::anyMissing(runValue(strand))) {
        warning("missing values in strand converted to \"*\"")
        runValue(strand)[is.na(runValue(strand))] <- "*"
    }

    lx <- max(length(seqnames), length(ranges), length(strand))
    if (lx > 1) {
        if (length(seqnames) == 1)
            seqnames <- rep(seqnames, lx)
        if (length(ranges) == 1)
            ranges <- rep(ranges, lx)
        if (length(strand) == 1)
            strand <- rep(strand, lx)
    }

    if (!is.integer(seqlengths))
        seqlengths <-
          structure(as.integer(seqlengths), names = names(seqlengths))

    elementMetadata <- DataFrame(...)
    if (ncol(elementMetadata) == 0)
        elementMetadata <- new("DataFrame", nrows = length(seqnames))
    if (!is.null(rownames(elementMetadata))) {
        if (!is.null(names(ranges)))
            names(ranges) <- rownames(elementMetadata)
        rownames(elementMetadata) <- NULL
    }

    new("GRanges", seqnames = seqnames, ranges = ranges, strand = strand,
        seqlengths = seqlengths, elementMetadata = elementMetadata)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("RangedData", "GRanges",
    function(from)
    {
        ranges <- unlist(ranges(from), use.names=FALSE)
        values <- unlist(values(from), use.names=FALSE)
        nms <- rownames(from)
        rownames(values) <- NULL
        whichStrand <- which(colnames(values) == "strand")
        if (length(whichStrand) > 0)
            values <- values[-whichStrand]
        GRanges(seqnames = space(from),
                ranges = ranges,
                strand = Rle(strand(from)),
                values)
    }
)

setAs("RangesList", "GRanges",
      function(from)
      {
        if (!length(from))
          return(GRanges())
        ranges <- unlist(from, use.names=FALSE)
        GRanges(seqnames = space(from),
                ranges = ranges,
                strand = Rle("*", length(ranges)),
                values = elementMetadata(ranges))
      })


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("seqnames", "GRanges", function(x) x@seqnames)
setMethod("ranges", "GRanges", function(x, ...) x@ranges)
setMethod("strand", "GRanges", function(x) x@strand)
setMethod("seqlengths", "GRanges", function(x) x@seqlengths)
setMethod("elementMetadata", "GRanges",
    function(x, ...)
    {
        ans <- x@elementMetadata
        if (!is.null(names(x)))
            rownames(ans) <- names(x)
        ans
    }
)

setReplaceMethod("seqnames", "GRanges",
    function(x, value) 
    {
        if (!is(value, "Rle"))
            value <- Rle(value)
        if (!is.factor(runValue(value)))
            runValue(value) <- factor(runValue(value))
        n <- length(x)
        k <- length(value)
        if (k != n) {
            if ((k == 0) || (k > n) || (n %% k != 0))
                stop(k, " elements in value to replace ", n, "elements")
            value <- rep(value, length.out = n)
        }
        seqlengths <- seqlengths(x)
        if (identical(runLength(seqnames(x)), runLength(value))) {
            matchTable <-
              unique(data.frame(old = runValue(seqnames(x)),
                                new = runValue(value)))
            if (!anyDuplicated(matchTable[["old"]]) &&
                !anyDuplicated(matchTable[["new"]])) {
                if (isTRUE(all.equal(as.integer(matchTable[["old"]]),
                                     as.integer(matchTable[["new"]])))) {
                    names(seqlengths) <- levels(value)
                } else {
                    names(seqlengths) <-
                        matchTable[["new"]][match(names(seqlengths),
                                            matchTable[["old"]])]
                }
            }
        }
        initialize(x, seqnames = value, seqlengths = seqlengths)
    }
)
setReplaceMethod("ranges", "GRanges",
    function(x, value) 
    {
        if (!is(value, "IRanges"))
            value <- as(value, "IRanges")
        n <- length(x)
        k <- length(value)
        if (k != n) {
            if ((k == 0) || (k > n) || (n %% k != 0))
                stop(k, " elements in value to replace ", n, "elements")
            value <- rep(value, length.out = n)
        }
        initialize(x, ranges = value)
    }
)
setReplaceMethod("strand", "GRanges",
    function(x, value) 
    {
        if (!is(value, "Rle"))
            value <- Rle(value)
        if (!is.factor(runValue(value)) ||
            !identical(levels(runValue(value)), levels(strand())))
            runValue(value) <- strand(runValue(value))
        n <- length(x)
        k <- length(value)
        if (k != n) {
            if ((k == 0) || (k > n) || (n %% k != 0))
                stop(k, " elements in value to replace ", n, "elements")
            value <- rep(value, length.out = n)
        }
        initialize(x, strand = value)
    }
)

setReplaceMethod("seqlengths", "GRanges",
    function(x, value)
    {
        if (!is.integer(value)) {
            nms <- names(value)
            value <- as.integer(value)
            names(value) <- nms
        }
        if (is.null(names(value))) {
            names(value) <- names(seqlengths(x))
        }
        if (length(setdiff(names(value), names(seqlengths(x)))) == 0) {
            initialize(x, seqlengths = value[names(seqlengths(x))])
        } else {
            seqnames <- seqnames(x)
            runValue(seqnames) <-
              factor(as.character(runValue(seqnames)), levels = names(value))
            initialize(x, seqnames = seqnames, seqlengths = value)
        }
    }
)

setReplaceMethod("elementMetadata", "GRanges",
    function(x, value)
    {
        if (is.null(value))
            value <- new("DataFrame", nrows = length(x))
        else if (!is(value, "DataFrame"))
            value <- DataFrame(value)
        if (!is.null(rownames(value)))
            rownames(value) <- NULL
        n <- length(x)
        k <- nrow(value)
        if (k != n) {
            if ((k == 0) || (k > n) || (n %% k != 0))
                stop(k, " rows in value to replace ", n, "rows")
            value <- value[rep(seq_len(k), length.out = n), , drop=FALSE]
        }
        initialize(x, elementMetadata = value)
    }
)

setMethod("names", "GRanges", function(x) names(x@ranges))
setReplaceMethod("names", "GRanges",
    function(x, value)
    {
        names(x@ranges) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ranges methods.
###

setReplaceMethod("start", "GRanges",
    function(x, check = TRUE, value)
    {
        if (!is.integer(value))
            value <- as.integer(value)
        ranges <- ranges(x)
        starts <- start(ranges)
        starts[] <- value
        if (!IRanges:::anyMissing(seqlengths(x))) {
            if (IRanges:::anyMissingOrOutside(starts, 1L)) {
                warning("trimmed start values to be positive")
                starts[starts < 1L] <- 1L
            }
        }
        start(ranges, check = check) <- starts
        initialize(x, ranges = ranges)
    }
)

setReplaceMethod("end", "GRanges",
    function(x, check = TRUE, value)
    {
        if (!is.integer(value))
            value <- as.integer(value)
        ranges <- ranges(x)
        ends <- end(ranges)
        ends[] <- value
        seqlengths <- seqlengths(x)
        if (!IRanges:::anyMissing(seqlengths)) {
            seqlengths <- seqlengths[levels(seqnames(x))]
            maxEnds <- seqlengths[as.integer(seqnames(x))]
            trim <- which(ends > maxEnds)
            if (length(trim) > 0) {
                warning("trimmed end values to be <= seqlengths")
                ends[trim] <- maxEnds[trim]
            }
        }
        end(ranges, check = check) <- ends
        initialize(x, ranges = ranges)
    }
)

setReplaceMethod("width", "GRanges",
    function(x, check = TRUE, value)
    {
        if (!is.integer(value))
            value <- as.integer(value)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            end(x) <- start(x) + (value - 1L)
        } else {
            ranges <- ranges(x)
            width(ranges, check = check) <- value
            x <- initialize(x, ranges = ranges)
        }
        x
    }
)

setMethod("flank", "GRanges",
    function(x, width, start = TRUE, both = FALSE, use.names = TRUE)
    {
        start <- as.vector(start == (strand(x) != "-"))
        ranges <-
          flank(ranges(x), width = width, start = start, both = both,
                use.names = use.names)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            start(x) <- start(ranges)
            end(x) <- end(ranges)
        } else {
            x <- initialize(x, ranges = ranges)
        }
        x
    }
)

setMethod("resize", "GRanges",
    function(x, width, fix = "start", use.names = TRUE)
    {
        if (!identical(fix, "start"))
            stop("'fix' arguments is not supported for GRanges objects")
        fix <- x@strand
        levels(fix) <- c("start", "end", "center")
        runValue(fix) <- as.character(runValue(fix))
        ranges <-
          resize(ranges(x), width = width, fix = fix, use.names = use.names)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            start(x) <- start(ranges)
            end(x) <- end(ranges)
        } else {
            x <- initialize(x, ranges = ranges)
        }
        x
    }
)

setMethod("shift", "GRanges",
    function(x, shift, use.names = TRUE)
    {
        ranges <- shift(ranges(x), shift, use.names = use.names)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            end(x) <- end(ranges)
            start(x) <- pmin.int(start(ranges), end(x))
        } else {
            x <- initialize(x, ranges = ranges)
        }
        x
    }
)

.interIntervalGRanges <- function(x, FUN, ...)
{
    xIRangesList <-
      split(unname(ranges(x)), paste(seqnames(x), strand(x), sep = "\r"))
    ansIRangesList <- FUN(xIRangesList, ...)
    k <- elementLengths(ansIRangesList)
    splitListNames <- strsplit(names(ansIRangesList), split = "\r")
    ansSeqnames <-
      Rle(factor(unlist(lapply(splitListNames, "[[", 1L)),
                 levels = levels(seqnames(x))), k)
    ansStrand <- Rle(strand(unlist(lapply(splitListNames, "[[", 2L))), k)
    GRanges(seqnames = ansSeqnames,
            ranges = unlist(ansIRangesList, use.names=FALSE),
            strand = ansStrand,
            seqlengths = seqlengths(x))
}

setMethod("disjoin", "GRanges",
    function(x)
        .interIntervalGRanges(x, disjoin)
)

setMethod("gaps", "GRanges",
    function(x, start = 1L, end = seqlengths(x))
    {
        seqlevels <- levels(seqnames(x))
        if (!is.null(names(start)))
            start <- start[seqlevels]
        if (!is.null(names(end)))
            end <- end[seqlevels]
        start <- IRanges:::recycleVector(start, length(seqlevels))
        start <- rep(start, each = 3)
        end <- IRanges:::recycleVector(end, length(seqlevels))
        end <- rep(end, each = 3)
        .interIntervalGRanges(x, gaps, start = start, end = end)
    }
)

setMethod("range", "GRanges",
    function(x, ..., na.rm)
        .interIntervalGRanges(unname(c(x, ...)), range)
)

setMethod("reduce", "GRanges",
    function(x, drop.empty.ranges = FALSE, min.gapwidth = 1L,
             with.inframe.attrib = FALSE)
    {
        if (!identical(with.inframe.attrib, FALSE))
            stop("'with.inframe.attrib' argument not supported ",
                 "when reducing a GRanges object")
        .interIntervalGRanges(x, reduce, drop.empty.ranges = drop.empty.ranges,
                              min.gapwidth = min.gapwidth)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sequence methods.
###

setMethod("[", "GRanges",
    function(x, i, j, ..., drop)
    {
        if (missing(i)) {
            if (!missing(j))
                x <-
                  initialize(x,
                             elementMetadata =
                             x@elementMetadata[, j, drop=FALSE])
        } else {
            iInfo <- IRanges:::.bracket.Index(i, length(x), names(x))
            if (!is.null(iInfo[["msg"]]))
                stop(iInfo[["msg"]])
            i <- iInfo[["idx"]]
            if (missing(j))
                elementMetadata <- x@elementMetadata[i, , drop=FALSE]
            else
                elementMetadata <- x@elementMetadata[i, j, drop=FALSE]
            ranges <- x@ranges[i]
            nms <- names(ranges)
            if (!is.null(nms)) {
                whichEmpty <- which(nms == "")
                nms[whichEmpty] <- as.character(whichEmpty)
                nms2 <- make.unique(nms)
                if (length(whichEmpty) > 0 || !identical(nms, nms2))
                    names(ranges) <- nms2
            }
            x <-
              initialize(x,
                         seqnames = x@seqnames[i],
                         ranges = ranges,
                         strand = x@strand[i],
                         elementMetadata = elementMetadata)
        }
        x
    }
)

setReplaceMethod("[", "GRanges",
    function(x, i, j, ..., value)
    {
        if (!is(value, "GRanges"))
            stop("replacement value must be a GRanges object")
        if (!identical(seqlengths(x), seqlengths(value)))
            stop("'seqlengths(x)' and 'seqlengths(value)' are not identical")
        seqnames <- x@seqnames
        ranges <- x@ranges
        strand <- x@strand
        elementMetadata <- x@elementMetadata
        if (missing(i)) {
            seqnames[] <- value@seqnames
            ranges[] <- value@ranges
            strand[] <- value@strand
            if (missing(j))
                elementMetadata[,] <- value@elementMetadata
            else
                elementMetadata[,j] <- value@elementMetadata
        } else {
            iInfo <- IRanges:::.bracket.Index(i, length(x), names(x))
            if (!is.null(iInfo[["msg"]]))
                stop(iInfo[["msg"]])
            i <- iInfo[["idx"]]
            seqnames[i] <- value@seqnames
            ranges[i] <- value@ranges
            strand[i] <- value@strand
            if (missing(j))
                elementMetadata[i,] <- value@elementMetadata
            else
                elementMetadata[i,j] <- value@elementMetadata
        }
        initialize(x, seqnames = seqnames, ranges = ranges,
                   strand = strand, elementMetadata = elementMetadata)
    }
)

setMethod("c", "GRanges",
    function(x, ..., recursive = FALSE)
    {
        if (recursive)
            stop("'recursive' mode not supported")
        args <- unname(list(x, ...))
        seqnames <- do.call(c, lapply(args, slot, "seqnames"))
        seqlengths <-
          do.call(c, lapply(args, slot, "seqlengths"))[levels(seqnames)]
        ranges <- do.call(c, lapply(args, slot, "ranges"))
        nms <- names(ranges)
        if (!is.null(nms)) {
            whichEmpty <- which(nms == "")
            nms[whichEmpty] <- as.character(whichEmpty)
            nms2 <- make.unique(nms)
            if (length(whichEmpty) > 0 || !identical(nms, nms2))
                names(ranges) <- nms2
        }
        initialize(x,
                   seqnames = seqnames,
                   ranges = ranges,
                   strand = do.call(c, lapply(args, slot, "strand")),
                   seqlengths = seqlengths,
                   elementMetadata =
                   do.call(rbind, lapply(args, slot, "elementMetadata")))
    }
)

setMethod("rev", "GRanges",
    function(x)
    {
        if (length(x) == 0)
            x
        else
            initialize(x, seqnames = rev(x@seqnames), ranges = rev(x@ranges),
                       strand = rev(x@strand),
                       elementMetadata =
                       x@elementMetadata[length(x):1, , drop=FALSE])
    }
)

setMethod("seqselect", "GRanges",
    function(x, start = NULL, end = NULL, width = NULL)
    {
        if (!is.null(end) || !is.null(width))
            start <- IRanges(start = start, end = end, width = width)
        irInfo <-
          IRanges:::.bracket.Index(start, length(x), names(x), asRanges = TRUE)
        if (!is.null(irInfo[["msg"]]))
            stop(irInfo[["msg"]])
        if (irInfo[["useIdx"]]) {
            ir <- irInfo[["idx"]]
            ranges <- seqselect(x@ranges, ir)
            nms <- names(ranges)
            if (!is.null(nms)) {
                whichEmpty <- which(nms == "")
                nms[whichEmpty] <- as.character(whichEmpty)
                nms2 <- make.unique(nms)
                if (length(whichEmpty) > 0 || !identical(nms, nms2))
                    names(ranges) <- nms2
            }
            x <-
              initialize(x,
                         seqnames = seqselect(x@seqnames, ir),
                         ranges = ranges,
                         strand = seqselect(x@strand, ir),
                         elementMetadata = seqselect(x@elementMetadata, ir))
        }
        x
    }
)

setReplaceMethod("seqselect", "GRanges",
    function(x, start = NULL, end = NULL, width = NULL, value)
    {
        if (!is(value, "GRanges"))
            stop("replacement value must be a GRanges object")
        if (!identical(seqlengths(x), seqlengths(value)))
            stop("'seqlengths(x)' and 'seqlengths(value)' are not identical")
        if (is.null(end) && is.null(width)) {
            if (is.null(start))
                ir <- IRanges(start = 1, width = length(x))
            else if (is(start, "Ranges"))
                ir <- start
            else {
                if (is.logical(start) && length(start) != length(x))
                    start <- rep(start, length.out = length(x))
                ir <- as(start, "IRanges")
            }
        } else {
            ir <- IRanges(start=start, end=end, width=width, names=NULL)
        }
        ir <- reduce(ir)
        if (length(ir) == 0) {
            x
        } else {
            seqnames <- as.vector(x@seqnames)
            ranges <- x@ranges
            strand <- as.vector(x@strand)
            elementMetadata <- x@elementMetadata
            seqselect(seqnames, ir) <- as.vector(value@seqnames)
            seqselect(ranges, ir) <- value@ranges
            seqselect(strand, ir) <- as.vector(value@strand)
            seqselect(elementMetadata, ir) <- value@elementMetadata
            initialize(x, seqnames = Rle(seqnames), ranges = ranges, 
                       strand = Rle(strand), elementMetadata = elementMetadata)
        }
    }
)

setMethod("split", "GRanges",
    function(x, f, drop = FALSE)
    {
        if (missing(f)) {
            nms <- names(x)
            if (is.null(nms))
                f <- seq_len(length(x))
            else
                f <- factor(nms, levels = nms)
        }
        nrows <- nlevels(f)
        if (nrows == 0)
            nrows <- sum(!is.na(unique(f)))
        IRanges:::newCompressedList("GRangesList", x, splitFactor = f,
                                    drop = drop,
                                    elementMetadata =
                                    new("DataFrame", nrows = nrows))
    }
)

setMethod("window", "GRanges",
    function(x, start = NA, end = NA, width = NA,
             frequency = NULL, delta = NULL, ...)
    {
        initialize(x,
                   seqnames =
                   window(x@seqnames, start = start, end = end, width = width,
                          frequency = frequency, delta = delta),
                   ranges =
                   window(x@ranges, start = start, end = end, width = width,
                          frequency = frequency, delta = delta),
                   strand =
                   window(x@strand, start = start, end = end, width = width,
                          frequency = frequency, delta = delta),
                   elementMetadata =
                   window(x@elementMetadata, start = start, end = end,
                          width = width, frequency = frequency, delta = delta))
    }
)

