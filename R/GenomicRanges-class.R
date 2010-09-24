### =========================================================================
### The GenomicRanges interface
### -------------------------------------------------------------------------
###

setClass("GenomicRanges", contains = c("Sequence", "VIRTUAL"))

### The code in this file will work out-of-the-box on 'x' as long as
### seqnames(x), ranges(x), strand(x), seqlengths(x), seqinfo(),
### update(x) and clone(x) are defined.

### TODO: a GenomicRangesList would be nice

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
        "'seqnames' should be a 'factor' Rle"
    else if (IRanges:::anyMissing(runValue(seqnames(x))))
        "'seqnames' contains missing values"
    else
        NULL
}

.valid.GenomicRanges.strand <- function(x)
{
    if (!is.factor(runValue(strand(x))) ||
        !identical(levels(runValue(strand(x))), levels(strand())))
        paste("'strand' should be a 'factor' Rle with levels c(",
              paste('"', levels(strand()), '"', sep = "", collapse = ", "),
                    ")", sep = "")
    else if (IRanges:::anyMissing(runValue(strand(x))))
        "'strand' contains missing values"
    else
        NULL
}

.valid.GenomicRanges.seqlengths <- function(x)
{
    msg <- NULL
    if (length(x) > 0 && is.null(names(seqlengths(x))))
        msg <- "'seqlengths' is unnamed"
    if (!setequal(names(seqlengths(x)), levels(seqnames(x))))
        msg <-
          c(msg, "'seqlengths' names do not match 'levels(seqnames)'")
    if (IRanges:::anyMissing(seqlengths(x))) {
        if (!all(is.na(seqlengths(x))))
            msg <-
              c(msg,
                "'seqlengths' cannot mix NAs and non-NA integers")
    } else if (setequal(names(seqlengths(x)), levels(seqnames(x)))) {
        if (any(seqlengths(x) < 0L, na.rm = TRUE))
            msg <- c(msg, "'seqlengths' contains negative values")
        ## TODO: Loosen the check below for circular sequences.
        seqnames <- seqnames(x)
        runValue(seqnames) <- runValue(seqnames)[drop=TRUE]
        minStarts <- IRanges:::.tapplyDefault(start(x), seqnames, min)
        maxEnds <- IRanges:::.tapplyDefault(end(x), seqnames, max)
        if (any(minStarts < 1L) || any(maxEnds > seqlengths(x)[names(maxEnds)]))
            msg <-
              c(msg, "'ranges' contains values outside of sequence bounds")
    }
    msg
}

.valid.GenomicRanges.isCircular <- function(x)
{
    msg <- NULL
    if (length(x) > 0 && is.null(names(isCircular(x))))
        msg <- "'isCircular' is unnamed"
    if (!setequal(names(isCircular(x)), levels(seqnames(x))))
        msg <-
          c(msg, "'isCircular' names do not match 'levels(seqnames)'")
    if (IRanges:::anyMissing(isCircular(x))) {
        if (!all(is.na(isCircular(x))))
            msg <-
              c(msg,
                "'isCircular' cannot mix NAs and non-NA integers")
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
    function(x, row.names = TRUE, ...)
    {
        if (!IRanges:::isTRUEorFALSE(row.names))
          stop("'row.names' must be TRUE or FALSE")
        ans <- callNextMethod()
        if (!row.names)
          rownames(ans) <- NULL
        ans
    }
)

setMethod("names", "GenomicRanges", function(x) names(ranges(x)))

setReplaceMethod("seqnames", "GenomicRanges",
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
                   ## Safe because 'names(seqlengths)' is guaranteed
                   ## to match the rows in 'seqinfo(x)'. Otherwise, we
                   ## would need to use Seqinfo().
                   seqinfo <- initialize(seqinfo(x),
                                         seqnames = names(seqlengths))
                   update(x, seqnames = value, seqinfo = seqinfo)
                 }
                 )
setReplaceMethod("ranges", "GenomicRanges",
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
                   update(x, ranges = value)
                 }
                 )
setReplaceMethod("strand", "GenomicRanges",
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
                   update(x, strand = value)
                 }
                 )

setReplaceMethod("seqlengths", "GenomicRanges",
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
                     seqlengths <- value[names(seqlengths(x))]
                     ## Safe because 'seqlengths' is guaranteed to
                     ## match the rows in 'seqinfo(x)'. Otherwise, we
                     ## would need to use Seqinfo().
                     seqinfo <- initialize(seqinfo(x), seqlengths = seqlengths)
                     update(x, seqinfo = seqinfo)
                   } else {
                     seqnames <- seqnames(x)
                     runValue(seqnames) <-
                       factor(as.character(runValue(seqnames)), levels = names(value))
                     is_circular <- isCircular(x)[names(value)]
                     ## The 'initialize(seqinfo(x), ...)' form would
                     ## not be safe here because we are resizing
                     ## 'seqinfo(x)'. Need to use Seqinfo() to recreate
                     ## the object from scratch.
                     seqinfo <- Seqinfo(seqnames = names(value),
                                        seqlengths = value,
                                        isCircular = is_circular)
                     update(x, seqnames = seqnames, seqinfo = seqinfo)
                   }
                 }
                 )

setReplaceMethod("elementMetadata", "GenomicRanges",
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
                   update(x, elementMetadata = value)
                 }
                 )

setReplaceMethod("names", "GenomicRanges",
                 function(x, value)
                 {
                   names(ranges(x)) <- value
                   x
                 }
                 )

setGeneric("seqinfo", function(x, ...) standardGeneric("seqinfo"))

setMethod("seqlengths", "GenomicRanges", function(x) seqlengths(seqinfo(x)))
setMethod("isCircular", "GenomicRanges", function(x) isCircular(seqinfo(x)))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning.
###
## An object is either 'update'd in place (usually with a replacement
## method) or 'clone'd (copied), with specified slots/fields overridden.

## For an object with a pure S4 slot representation, these both map to
## initialize. Reference classes will want to override 'update'. Other
## external representations need further customization.

setMethod("update", "ANY",
          function(object, ...)
          {
            initialize(object, ...)
          })

setGeneric("clone", function(x, ...) standardGeneric("clone"))
setMethod("clone", "ANY",
          function(x, ...)
          {
            initialize(x, ...)
          })

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
                                          shift = shift[[i]],
                                          width = width[[i]],
                                          weight = weight[[i]])
                             }))
    }
)

setReplaceMethod("start", "GenomicRanges",
                 function(x, check = TRUE, value)
                 {
                   if (!is.integer(value))
                     value <- as.integer(value)
                   ranges <- ranges(x)
                   starts <- start(ranges)
                   starts[] <- value
                   ## TODO: Revisit this to handle circularity (maybe).
                   if (!IRanges:::anyMissing(seqlengths(x))) {
                     if (IRanges:::anyMissingOrOutside(starts, 1L)) {
                       warning("trimmed start values to be positive")
                       starts[starts < 1L] <- 1L
                     }
                   }
                   start(ranges, check = check) <- starts
                   update(x, ranges = ranges)
                 }
                 )

setReplaceMethod("end", "GenomicRanges",
                 function(x, check = TRUE, value)
                 {
                   if (!is.integer(value))
                     value <- as.integer(value)
                   ranges <- ranges(x)
                   ends <- end(ranges)
                   ends[] <- value
                   seqlengths <- seqlengths(x)
                   ## TODO: Revisit this to handle circularity.
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
                   update(x, ranges = ranges)
                 }
                 )

setReplaceMethod("width", "GenomicRanges",
                 function(x, check = TRUE, value)
                 {
                   if (!is.integer(value))
                     value <- as.integer(value)
                   if (!IRanges:::anyMissing(seqlengths(x))) {
                     end(x) <- start(x) + (value - 1L)
                   } else {
                     ranges <- ranges(x)
                     width(ranges, check = check) <- value
                     x <- update(x, ranges = ranges)
                   }
                   x
                 }
                 )

setMethod("flank", "GenomicRanges",
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
              x <- clone(x, ranges = ranges)
            }
            x
          }
          )

setMethod("resize", "GenomicRanges",
          function(x, width, fix = "start", use.names = TRUE)
          {
            if (!missing(fix) &&
                (length(fix) > length(x) || length(x) %% length(fix) > 0))
              stop("'x' is not a multiple of 'fix' length")
            revFix <- c(start = "end", end = "start", center = "center")
            fix <- ifelse(strand(x) == "-", revFix[fix], fix)
            ranges <-
              resize(ranges(x), width = width, fix = fix, use.names = use.names)
            if (!IRanges:::anyMissing(seqlengths(x))) {
              start(x) <- start(ranges)
              end(x) <- end(ranges)
            } else {
              x <- clone(x, ranges = ranges)
            }
            x
          }
          )

setMethod("shift", "GenomicRanges",
          function(x, shift, use.names = TRUE)
          {
            ranges <- shift(ranges(x), shift, use.names = use.names)
            if (!IRanges:::anyMissing(seqlengths(x))) {
              end(x) <- end(ranges)
              start(x) <- pmin.int(start(ranges), end(x))
            } else {
              x <- clone(x, ranges = ranges)
            }
            x
          }
          )

.interIntervalGenomicRanges <- function(x, FUN, ...)
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
  clone(x, seqnames = ansSeqnames,
        ranges = unlist(ansIRangesList, use.names=FALSE),
        strand = ansStrand,
        elementMetadata = new("DataFrame", nrows = length(ansSeqnames)))
}

setMethod("disjoin", "GenomicRanges",
          function(x)
          .interIntervalGenomicRanges(x, disjoin)
          )

setMethod("gaps", "GenomicRanges",
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
            .interIntervalGenomicRanges(x, gaps, start = start, end = end)
          }
          )

setMethod("range", "GenomicRanges",
          function(x, ..., na.rm)
          .interIntervalGenomicRanges(unname(c(x, ...)), range)
          )

setMethod("reduce", "GenomicRanges",
          function(x, drop.empty.ranges = FALSE, min.gapwidth = 1L,
                   with.inframe.attrib = FALSE)
          {
            if (!identical(with.inframe.attrib, FALSE))
              stop("'with.inframe.attrib' argument not supported ",
                   "when reducing a GenomicRanges object")
            .interIntervalGenomicRanges(x, reduce,
                                        drop.empty.ranges = drop.empty.ranges,
                                        min.gapwidth = min.gapwidth)
          }
          )

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sequence methods.
###

setMethod("length", "GenomicRanges", function(x) length(seqnames(x)))

setMethod("[", "GenomicRanges",
          function(x, i, j, ..., drop)
          {
            if (missing(i)) {
              if (!missing(j))
                x <-
                  clone(x,
                         elementMetadata =
                         elementMetadata(x, FALSE)[, j, drop=FALSE])
            } else {
              iInfo <- IRanges:::.bracket.Index(i, length(x), names(x))
              if (!is.null(iInfo[["msg"]]))
                stop(iInfo[["msg"]])
              i <- iInfo[["idx"]]
              if (missing(j))
                elementMetadata <- elementMetadata(x, FALSE)[i, , drop=FALSE]
              else
                elementMetadata <- elementMetadata(x, FALSE)[i, j, drop=FALSE]
              ranges <- ranges(x)[i]
              nms <- names(ranges)
              if (!is.null(nms)) {
                whichEmpty <- which(nms == "")
                nms[whichEmpty] <- as.character(whichEmpty)
                nms2 <- make.unique(nms)
                if (length(whichEmpty) > 0 || !identical(nms, nms2))
                  names(ranges) <- nms2
              }
              x <-
                clone(x,
                      seqnames = seqnames(x)[i],
                      ranges = ranges,
                      strand = strand(x)[i],
                      elementMetadata = elementMetadata)
            }
            x
          }
          )

setReplaceMethod("[", "GenomicRanges",
                 function(x, i, j, ..., value)
                 {
                   if (!is(value, "GenomicRanges"))
                     stop("replacement value must be a GenomicRanges object")
                   ## TODO: Shouldn't we also compare the circularity flags?
                   if (!identical(seqlengths(x), seqlengths(value)))
                     stop("'seqlengths(x)' and 'seqlengths(value)' are not identical")
                   seqnames <- seqnames(x)
                   ranges <- ranges(x)
                   strand <- strand(x)
                   elementMetadata <- elementMetadata(x, FALSE)
                   if (missing(i)) {
                     seqnames[] <- seqnames(value)
                     ranges[] <- ranges(value)
                     strand[] <- strand(value)
                     if (missing(j))
                       elementMetadata[,] <- elementMetadata(value, FALSE)
                     else
                       elementMetadata[,j] <- elementMetadata(value, FALSE)
                   } else {
                     iInfo <- IRanges:::.bracket.Index(i, length(x), names(x))
                     if (!is.null(iInfo[["msg"]]))
                       stop(iInfo[["msg"]])
                     i <- iInfo[["idx"]]
                     seqnames[i] <- seqnames(value)
                     ranges[i] <- ranges(value)
                     strand[i] <- strand(value)
                     if (missing(j))
                       elementMetadata[i,] <- elementMetadata(value, FALSE)
                     else
                       elementMetadata[i,j] <- elementMetadata(value, FALSE)
                   }
                   update(x, seqnames = seqnames, ranges = ranges,
                          strand = strand, elementMetadata = elementMetadata)
                 }
                 )

setMethod("c", "GenomicRanges",
          function(x, ..., recursive = FALSE)
          {
            if (recursive)
              stop("'recursive' mode not supported")
            args <- unname(list(x, ...))
            seqnames <- do.call(c, lapply(args, seqnames))
            ## TODO: Revisit the 2 lines below. This code is silently ignoring
            ## the fact that we might be combining objects with incompatible
            ## sequence lengths and/or circularity flags. Is it reasonable?
            ## Note that it is inconsistent with what [<- does.
            seqlengths <- do.call(c, lapply(args, seqlengths))[levels(seqnames)]
            is_circular <-
              do.call(c, lapply(args, isCircular))[levels(seqnames)]
            ## The 'initialize(seqinfo(x), ...)' form would not be safe here
            ## because we are resizing 'seqinfo(x)'. Need to use Seqinfo() to
            ## recreate the object from scratch.
            seqinfo <- Seqinfo(seqnames = names(seqlengths),
                               seqlengths = seqlengths,
                               isCircular = is_circular)
            ranges <- do.call(c, lapply(args, ranges))
            nms <- names(ranges)
            if (!is.null(nms)) {
              whichEmpty <- which(nms == "")
              nms[whichEmpty] <- as.character(whichEmpty)
              nms2 <- make.unique(nms)
              if (length(whichEmpty) > 0 || !identical(nms, nms2))
                names(ranges) <- nms2
            }
            clone(x,
                  seqnames = seqnames,
                  ranges = ranges,
                  strand = do.call(c, lapply(args, strand)),
                  seqinfo = seqinfo,
                  elementMetadata =
                  do.call(rbind, lapply(args, elementMetadata, FALSE)))
          }
          )

setMethod("rev", "GenomicRanges",
          function(x)
          {
            if (length(x) == 0)
              x
            else
              clone(x, seqnames = rev(seqnames(x)),
                    ranges = rev(ranges(x)),
                    strand = rev(strand(x)),
                    elementMetadata =
                    elementMetadata(x, FALSE)[length(x):1, , drop=FALSE])
          }
          )

setMethod("seqselect", "GenomicRanges",
          function(x, start = NULL, end = NULL, width = NULL)
          {
            if (!is.null(end) || !is.null(width))
              start <- IRanges(start = start, end = end, width = width)
            irInfo <-
              IRanges:::.bracket.Index(start, length(x), names(x),
                                       asRanges = TRUE)
            if (!is.null(irInfo[["msg"]]))
              stop(irInfo[["msg"]])
            if (irInfo[["useIdx"]]) {
              ir <- irInfo[["idx"]]
              ranges <- seqselect(ranges(x), ir)
              nms <- names(ranges)
              if (!is.null(nms)) {
                whichEmpty <- which(nms == "")
                nms[whichEmpty] <- as.character(whichEmpty)
                nms2 <- make.unique(nms)
                if (length(whichEmpty) > 0 || !identical(nms, nms2))
                  names(ranges) <- nms2
              }
              x <-
                clone(x,
                      seqnames = seqselect(seqnames(x), ir),
                      ranges = ranges,
                      strand = seqselect(strand(x), ir),
                      elementMetadata =
                      seqselect(elementMetadata(x, FALSE), ir))
            }
            x
          }
          )

setReplaceMethod("seqselect", "GenomicRanges",
                 function(x, start = NULL, end = NULL, width = NULL, value)
                 {
                   if (!is(value, "GenomicRanges"))
                     stop("replacement value must be a GenomicRanges object")
                   if (!identical(seqlengths(x), seqlengths(value)))
                     stop("'seqlengths(x)' and 'seqlengths(value)' are not identical")
                   ## TODO: Shouldn't we also compare the circularity flags?
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
                     ir <- IRanges(start=start, end=end, width=width,
                                   names=NULL)
                   }
                   ir <- reduce(ir)
                   if (length(ir) == 0) {
                     x
                   } else {
                     seqnames <- as.factor(seqnames(x))
                     ranges <- ranges(x)
                     strand <- as.factor(strand(x))
                     elementMetadata <- elementMetadata(x, FALSE)
                     seqselect(seqnames, ir) <- as.factor(seqnames(value))
                     seqselect(ranges, ir) <- ranges(value)
                     seqselect(strand, ir) <- as.factor(strand(value))
                     seqselect(elementMetadata, ir) <- elementMetadata(value)
                     update(x, seqnames = Rle(seqnames), ranges = ranges, 
                            strand = Rle(strand),
                            elementMetadata = elementMetadata)
                   }
                 }
                 )

setMethod("window", "GenomicRanges",
          function(x, start = NA, end = NA, width = NA,
                   frequency = NULL, delta = NULL, ...)
          {
            update(x,
                   seqnames =
                   window(seqnames(x), start = start, end = end,
                          width = width, frequency = frequency,
                          delta = delta),
                   ranges =
                   window(ranges(x), start = start, end = end,
                          width = width, frequency = frequency,
                          delta = delta),
                   strand =
                   window(strand(x), start = start, end = end,
                          width = width, frequency = frequency,
                          delta = delta),
                   elementMetadata =
                   window(elementMetadata(x, FALSE), start = start, end = end,
                          width = width, frequency = frequency,
                          delta = delta))
          }
          )


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

