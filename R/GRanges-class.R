### =========================================================================
### GRanges objects
### -------------------------------------------------------------------------
###
### Class definition

setClass("GRanges", contains = "Sequence",
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

.valid.GRanges.slots <- function(x)
{
    n <- length(x@seqnames)
    if ((length(x@ranges) != n) || (length(x@strand) != n) ||
        (nrow(x@elementMetadata) != n))
        "slot lengths are not all equal"
    else
        NULL
}

.valid.GRanges.seqnames <- function(x)
{
    if (!is.factor(runValue(x@seqnames)))
        "slot 'seqnames' should be a 'factor' Rle"
    else if (IRanges:::anyMissing(runValue(x@seqnames)))
        "slot 'seqnames' contains missing values"
    else
        NULL
}

.valid.GRanges.strand <- function(x)
{
    if (!is.factor(runValue(x@strand)) ||
        !identical(levels(runValue(x@strand)), levels(strand())))
        paste("slot 'strand' should be a 'factor' Rle with levels c(",
              paste('"', levels(strand()), '"', sep = "", collapse = ", "),
                    ")", sep = "")
    else if (IRanges:::anyMissing(runValue(x@strand)))
        "slot 'strand' contains missing values"
    else
        NULL
}

.valid.GRanges.seqlengths <- function(x)
{
    msg <- NULL
    if (length(x) > 0 && is.null(names(x@seqlengths)))
        msg <- "slot 'seqlengths' is unnamed"
    if (!setequal(names(x@seqlengths), levels(x@seqnames)))
        msg <-
          c(msg, "slot 'seqlengths' names to not match 'levels(seqnames)'")
    if (any(x@seqlengths < 0L, na.rm = TRUE))
        msg <- c(msg, "slot 'seqlengths' contains negative values")
    if (IRanges:::anyMissing(x@seqlengths)) {
        if (!all(is.na(x@seqlengths)))
            msg <-
              c(msg,
                "slot 'seqlengths' cannot mix NAs and non-NA integers")
    } else {
        seqnames <- seqnames(x)
        runValue(seqnames) <- runValue(seqnames)[drop=TRUE]
        minStarts <- IRanges:::.tapplyDefault(start(x), seqnames, min)
        maxEnds <- IRanges:::.tapplyDefault(end(x), seqnames, max)
        if (any(minStarts < 1L) || any(maxEnds > x@seqlengths[names(maxEnds)]))
            msg <-
              c(msg, "slot 'ranges' contains values outside of sequence bounds")
    }
    msg
}

.valid.GRanges.elementMetadata <- function(x)
{
    msg <- NULL
    if (any(c("seqnames", "ranges", "strand", "seqlengths", "start", "end",
              "width", "element") %in% colnames(x@elementMetadata)))
        msg <-
          paste("slot 'elementMetadata' cannot use \"seqnames\", \"ranges\",",
                "\"strand\", \"seqlengths\", \"start\", \"end\", \"width\", or",
                "\"element\" as column names")
    if (!is.null(rownames(x@elementMetadata)))
        msg <- c(msg, "slot 'elementMetadata' cannot contain row names")
    msg
}

.valid.GRanges <- function(x)
{
    c(.valid.GRanges.slots(x),
      .valid.GRanges.seqnames(x),
      .valid.GRanges.strand(x),
      .valid.GRanges.seqlengths(x),
      .valid.GRanges.elementMetadata(x))
}

setValidity2("GRanges", .valid.GRanges)


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

setMethod("as.data.frame", "GRanges",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        ranges <- ranges(x)
        if (missing(row.names))
            row.names <- names(x)
        if (!is.null(names(x)))
            names(x) <- NULL
        data.frame(seqnames = as.vector(seqnames(x)),
                   start = start(x),
                   end = end(x),
                   width = width(x),
                   strand = as.vector(strand(x)),
                   as.data.frame(elementMetadata(x)),
                   row.names = row.names,
                   stringsAsFactors = FALSE)
    }
)


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
            matchTable <- table(runValue(seqnames(x)), runValue(value)) > 0L
            if ((nrow(matchTable) == ncol(matchTable)) &&
                all(rowSums(matchTable) == 1L) &&
                all(colSums(matchTable) == 1L)) {
                names(seqlengths) <- colnames(matchTable)[max.col(matchTable)]
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
        if (!is.integer(value))
            value <- as.integer(value)
        initialize(x, seqlengths = value)
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

setMethod("start", "GRanges", function(x, ...) start(x@ranges))
setMethod("end", "GRanges", function(x, ...) end(x@ranges))
setMethod("width", "GRanges", function(x) width(x@ranges))

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
    function(x, start=NA, end=NA)
        .interIntervalGRanges(x, gaps, start = start, end = end)
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

setMethod("coverage", "GRanges",
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

setMethod("length", "GRanges", function(x) length(x@seqnames))

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
        IRanges:::newCompressedList("GRangesList", x, splitFactor = f,
                                    drop = drop,
                                    elementMetadata =
                                    new("DataFrame",
                                        nrows = sum(!is.na(unique(f)))))
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

.showGRanges <-
function(object, print.seqlengths = FALSE)
{
    lo <- length(object)
    nc <- ncol(object@elementMetadata)
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
}

setMethod("show", "GRanges",
          function(object) .showGRanges(object, print.seqlengths = TRUE))