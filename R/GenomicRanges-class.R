### =========================================================================
### The GenomicRanges interface
### -------------------------------------------------------------------------
###

setClass("GenomicRanges",
    contains="Sequence",
    representation("VIRTUAL")
)

### The code in this file will work out-of-the-box on 'x' as long as
### seqnames(x), ranges(x), strand(x), seqlengths(x), seqinfo(),
### update(x) and clone(x) are defined.

### TODO: a GenomicRangesList would be nice


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### 2 non-exported low-level helper functions.
###
### For the 2 functions below, 'x' must be a GenomicRanges object.
### They both return a named integer vector where the names are guaranteed
### to be 'seqlevels(x)'.
###

minStartPerGRangesSequence <- function(x)
{
    cil <- IRanges:::newCompressedList("CompressedIntegerList",
                       unlistData=start(x),
                       splitFactor=seqnames(x))
    v <- Views(cil@unlistData, cil@partitioning)  # XIntegerViews object
    ans <- viewMins(v)
    ans[width(v) == 0L] <- NA_integer_
    names(ans) <- names(v)
    ans
}

maxEndPerGRangesSequence <- function(x)
{
    cil <- IRanges:::newCompressedList("CompressedIntegerList",
                       unlistData=end(x),
                       splitFactor=seqnames(x))
    v <- Views(cil@unlistData, cil@partitioning)  # XIntegerViews object
    ans <- viewMaxs(v)
    ans[width(v) == 0L] <- NA_integer_
    names(ans) <- names(v)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setMethod("length", "GenomicRanges", function(x) length(seqnames(x)))

setMethod("names", "GenomicRanges", function(x) names(ranges(x)))

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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###
### TODO: Should we enforce ranges(x) to be an IRanges *instance* (i.e.
### class(ranges(x) == "IRanges")) instead of just an IRanges *object* (i.e.
### is(ranges(x), "IRanges"))? Right now I can create a GRanges object where
### the ranges are a Views object on a very long DNAString subject with
### something like 'GRanges("chr1", Views(subject, start=1:2, end=5))'.
### Sounds cool but there are also some potential complications with this...

.valid.GenomicRanges.length <- function(x)
{
    n <- length(seqnames(x))
    if ((length(ranges(x)) != n)
     || (length(strand(x)) != n)
     || (nrow(elementMetadata(x)) != n))
        return("slot lengths are not all equal")
    NULL
}

.valid.GenomicRanges.seqnames <- function(x)
{
    if (!is.factor(runValue(seqnames(x))))
        return("'seqnames' should be a 'factor' Rle")
    if (IRanges:::anyMissing(runValue(seqnames(x))))
        return("'seqnames' contains missing values")
    NULL
}

.valid.GenomicRanges.strand <- function(x)
{
    if (!is.factor(runValue(strand(x))) ||
        !identical(levels(runValue(strand(x))), levels(strand())))
    {
        msg <- c("'strand' should be a 'factor' Rle with levels c(",
                 paste('"', levels(strand()), '"', sep = "", collapse = ", "),
                 ")")
        return(paste(msg, collapse=""))
    }
    if (IRanges:::anyMissing(runValue(strand(x))))
        return("'strand' contains missing values")
    NULL
}

.valid.GenomicRanges.elementMetadata <- function(x)
{
    ## NOTE: This list is also included in the man page for GRanges objects.
    ## Keep the 2 lists in sync!
    INVALID.COLNAMES <- c("seqnames", "ranges", "strand",
                          "seqlevels", "seqlengths", "isCircular",
                          "start", "end", "width", "element")
    if (any(INVALID.COLNAMES %in% colnames(elementMetadata(x)))) {
        msg <- c("slot 'elementMetadata' cannot use",
                 paste("\"", INVALID.COLNAMES, "\"", sep="", collapse=", "),
                 "as column names")
        return(paste(msg, collapse=" "))
    }
    NULL
}

### Also used by the validity method for GappedAlignments objects.
valid.GenomicRanges.seqinfo <- function(x)
{
    x_seqinfo <- seqinfo(x)
    if (!identical(seqlevels(x_seqinfo), levels(seqnames(x)))) {
        msg <- c("'seqlevels(seqinfo(x))' and 'levels(seqnames(x))'",
                 "are not identical")
        return(paste(msg, collapse=" "))
    }
    x_seqlengths <- seqlengths(x_seqinfo)
    seqs_with_known_length <- names(x_seqlengths)[!is.na(x_seqlengths)]
    if (length(seqs_with_known_length) == 0L)
        return(NULL)
    if (any(x_seqlengths < 0L, na.rm = TRUE))
        return("'seqlengths(x)' contains negative values")
    ## We check only the ranges that are on a non-circular sequence with
    ## a known length.
    x_isCircular <- isCircular(x_seqinfo)
    non_circ_seqs <- names(x_isCircular)[!(x_isCircular %in% TRUE)]
    ncswkl <- intersect(non_circ_seqs, seqs_with_known_length)
    if (length(ncswkl) == 0L)
        return(NULL)
    x_seqnames <- seqnames(x)
    runValue(x_seqnames) <- runValue(x_seqnames)[drop=TRUE]
    minStarts <- minStartPerGRangesSequence(x)
    maxEnds <- maxEndPerGRangesSequence(x)
    if (any(minStarts[ncswkl] < 1L, na.rm = TRUE)
     || any(maxEnds[ncswkl] >
            seqlengths(x)[ncswkl], na.rm = TRUE))
        return("'ranges' contains values outside of sequence bounds")
    NULL
}

.valid.GenomicRanges <- function(x)
{
    c(.valid.GenomicRanges.length(x),
      .valid.GenomicRanges.seqnames(x),
      .valid.GenomicRanges.strand(x),
      .valid.GenomicRanges.elementMetadata(x),
      valid.GenomicRanges.seqinfo(x))
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
        metadata(ranges(rd)) <- metadata(from)
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
        metadata(rl) <- metadata(from)
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
### Setters.
###

setReplaceMethod("names", "GenomicRanges",
    function(x, value)
    {
        names(ranges(x)) <- value
        x
    }
)

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
                                     as.integer(matchTable[["new"]]))))
                {
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

normargGenomicRangesStrand <- function(strand, n)
{
    if (!is(strand, "Rle"))
        strand <- Rle(strand)
    if (!is.factor(runValue(strand))
     || !identical(levels(runValue(strand)), levels(strand())))
        runValue(strand) <- strand(runValue(strand))
    k <- length(strand)
    if (k != n) {
        if (k != 1L && (k == 0L || k > n || n %% k != 0L))
            stop("supplied 'strand' has ", k, " elements (", n, " expected)")
        strand <- rep(strand, length.out = n)
    }
    strand
}

setReplaceMethod("strand", "GenomicRanges",
    function(x, value) 
    {
        value <- normargGenomicRangesStrand(value, length(x))
        update(x, strand = value)
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
            if ((k == 0L) || (k > n) || (n %% k != 0L))
                stop(k, " rows in value to replace ", n, "rows")
            value <- value[rep(seq_len(k), length.out = n), , drop=FALSE]
        }
        update(x, elementMetadata = value)
    }
)

setReplaceMethod("seqlevels", "GenomicRanges",
    function(x, value)
    {
        if (!is.character(value) || any(is.na(value)))
            stop("supplied 'seqlevels' must be a character vector with no NAs")
        x_seqnames <- seqnames(x)
        dangling_seqlevels <- setdiff(unique(x_seqnames), value)
        if (length(dangling_seqlevels))
            stop("Supplied 'seqlevels' would exclude data for ",
                 paste(dangling_seqlevels, collapse = ", "),
                 ". Consider subsetting first.")
        runValue(x_seqnames) <-
          factor(as.character(runValue(x_seqnames)), levels = value)
        seqlengths <- seqlengths(x)[value]
        is_circular <- isCircular(x)[value]
        ## The 'initialize(seqinfo(x), ...)' form would
        ## not be safe here because we are resizing
        ## 'seqinfo(x)'. Need to use Seqinfo() to recreate
        ## the object from scratch.
        seqinfo <- Seqinfo(seqnames = value,
                           seqlengths = unname(seqlengths),
                           isCircular = unname(is_circular))
        update(x, seqnames = x_seqnames, seqinfo = seqinfo)
    }
)

### Needed only for BioC 2.8 to redirect users that are trying to modify the
### levels of 'seqnames(x)' via this method.
### TODO: Remove in BioC 2.9.
setReplaceMethod("seqlengths", "GenomicRanges",
    function(x, value)
    {
        if (!is.null(names(value)) && !setequal(names(value), seqlevels(x)))
            stop("The names of the supplied 'seqlengths' don't match ",
                 "'seqlevels(x)' (aka the\n  sequence levels). ",
                 "Please use the \"seqlevels\" setter if your intention ",
                 "was to\n  modify the sequence levels (i.e. use something ",
                 "like 'seqlevels(x) <- value'\n  where 'value' is ",
                 "a character vector containing the new levels).")
        callNextMethod()
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning.
###
### An object is either 'update'd in place (usually with a replacement
### method) or 'clone'd (copied), with specified slots/fields overridden.

### For an object with a pure S4 slot representation, these both map to
### initialize. Reference classes will want to override 'update'. Other
### external representations need further customization.

setMethod("update", "GenomicRanges",
          function(object, ...)
          {
            initialize(object, ...)
          })

setGeneric("clone", function(x, ...) standardGeneric("clone"))
setMethod("clone", "GenomicRanges",
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

### 'circle.length' must be NA (if the underlying sequence is linear) or the
### length of the underlying circular sequence (integer vector of length 1
### with the name of the sequence).
### 'rg' must be an IRanges object.
.coverage.circle <- function(circle.length, rg, shift, width, weight)
{
    if (is.na(circle.length))
        return(coverage(rg, shift = shift, width = width, weight = weight))
    cvg <- fold(coverage(rg, weight = weight),
                circle.length, from = 1L - shift)
    if (is.null(width))
        return(cvg)
    if (width > length(cvg))
        stop("invalid width (", width, ") ",
             "for circular sequence ", names(circle.length))
    cvg[seq_len(width)]
}

setMethod("coverage", "GenomicRanges",
    function(x, shift = list(0L), width = as.list(seqlengths(x)),
             weight = list(1L))
    {
        if (any(start(x) < 1L))
            stop("'x' contains ranges starting before position 1. ",
                 "coverage() currently doesn't support this.")
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
                         rg <- seqselect(xRanges, xSplitRanges[[i]])
                         if (isCircular(x)[i] %in% TRUE)
                             circle.length <- seqlengths(x)[i]
                         else
                             circle.length <- NA
                         .coverage.circle(circle.length, rg,
                                          shift[[i]], width[[i]], weight[[i]])
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
              end(x, check=FALSE) <- end(ranges)
              start(x) <- pmin.int(start(ranges), end(x))
            } else {
              x <- clone(x, ranges = ranges)
            }
            x
          }
          )

.interIntervalGenomicRanges <- function(x, FUN, ...)
{
    elementMetadata(x) <- NULL
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

setMethod("[", "GenomicRanges",
          function(x, i, j, ..., drop)
          {
            if (!missing(i)) {
              iInfo <- IRanges:::.bracket.Index(i, length(x), names(x))
              if (!is.null(iInfo[["msg"]]))
                stop(iInfo[["msg"]])
            }
            if (missing(i) || !iInfo[["useIdx"]]) {
              if (!missing(j))
                x <-
                  clone(x,
                         elementMetadata =
                         elementMetadata(x, FALSE)[, j, drop=FALSE])
            } else {
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
          })

setReplaceMethod("[", "GenomicRanges",
    function(x, i, j, ..., value)
    {
        if (!is(value, "GenomicRanges"))
            stop("replacement value must be a GenomicRanges object")
        seqinfo(x) <- merge(seqinfo(x), seqinfo(value))
        seqnames <- seqnames(x)
        ranges <- ranges(x)
        strand <- strand(x)
        elementMetadata <- elementMetadata(x, FALSE)
        if (!missing(i)) {
            iInfo <- IRanges:::.bracket.Index(i, length(x), names(x))
            if (!is.null(iInfo[["msg"]]))
                stop(iInfo[["msg"]])
        }
        if (missing(i) || !iInfo[["useIdx"]]) {
            seqnames[] <- seqnames(value)
            ranges[] <- ranges(value)
            strand[] <- strand(value)
            if (missing(j))
                elementMetadata[,] <- elementMetadata(value, FALSE)
            else
                elementMetadata[,j] <- elementMetadata(value, FALSE)
        } else {
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
        ans_seqinfo <- do.call(merge, lapply(args, seqinfo))
        ans_seqnames <- do.call(c, lapply(args, seqnames))
        ans_ranges <- do.call(c, lapply(args, ranges))
        ans_strand <- do.call(c, lapply(args, strand))
        ans_elementMetadata <- do.call(rbind, lapply(args, elementMetadata, FALSE))
        ans_names <- names(ans_ranges)
        if (!is.null(ans_names)) {
            whichEmpty <- which(ans_names == "")
            ans_names[whichEmpty] <- as.character(whichEmpty)
            ans_names2 <- make.unique(ans_names)
            if (length(whichEmpty) > 0 || !identical(ans_names, ans_names2))
                names(ans_ranges) <- ans_names2
        }
        clone(x,
              seqnames = ans_seqnames,
              ranges = ans_ranges,
              strand = ans_strand,
              seqinfo = ans_seqinfo,
              elementMetadata = ans_elementMetadata)
    }
)

setMethod("seqselect", "GenomicRanges",
    function(x, start = NULL, end = NULL, width = NULL)
    {
        if (!is.null(end) || !is.null(width))
            start <- IRanges(start = start, end = end, width = width)
        irInfo <- IRanges:::.bracket.Index(start, length(x), names(x),
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
                if (length(whichEmpty) > 0L || !identical(nms, nms2))
                    names(ranges) <- nms2
            }
            x <- clone(x,
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
        seqinfo(x) <- merge(seqinfo(x), seqinfo(value))
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
            ir <- IRanges(start = start, end = end, width = width, names = NULL)
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

setMethod("precede", c("GenomicRanges", "GenomicRanges"),
  function(x, subject, ...)
  {
      if(all(as.character(strand(x)) == "*") && all(as.character(strand(subject)) == "*"))
        strand(x) <- strand(subject) <- "+"
      
      elementMetadata(x) <- list("posIndx" = seq_len(length(x)))
      elementMetadata(subject) <- list("posIndx" = seq_len(length(subject)))
      xLst <- split(x, seqnames(x), drop = FALSE)
      xLst <- xLst[sapply(xLst, length) > 0]
      xSeqNames <- names(xLst)
      subjectLst <- split(subject, seqnames(subject), drop = FALSE)[xSeqNames]  
      matchPos <- rep.int(NA_integer_, length(x))

      for (sq in xSeqNames) {
          x1Split <- split(xLst[[sq]], strand(xLst[[sq]]))
          x1Split <- x1Split[lapply(x1Split, length) > 0]
          s1 <- subjectLst[[sq]]
          for (st in names(x1Split)) {
              if( st == "+" ){
                  subSplit <- unlist(split(s1, strand(s1))[c("+", "*")])
                  ## call precede
                  res <- precede(ranges(x1Split[[st]]), ranges(subSplit))
                  indx <- !is.na(res)
                  res <- res[indx]
                  matchPos[elementMetadata(x1Split[[st]])$posIndx[indx]] <-
                  elementMetadata(subSplit)$posIndx[res]

              }else if(st == "-"){
                  subSplit <- unlist(split(s1, strand(s1))[c("-", "*")])
                  ## call follow
                  res <- follow(ranges(x1Split[[st]]), ranges(subSplit))
                  indx <- !is.na(res)
                  res <- res[indx]
                  matchPos[elementMetadata(x1Split[[st]])$posIndx[indx]] <-
                  elementMetadata(subSplit)$posIndx[res]
              }else if(st == "*"){
                  subSplit1 <- unlist(split(s1, strand(s1))[c("+")])
                  ## call precede
                  res1 <- precede(ranges(x1Split[[st]]), ranges(subSplit1))
                  indx1 <- !is.na(res1)
                  res1 <- res1[indx1]
                  k1 <- elementMetadata(x1Split[[st]])$posIndx[indx1]
                  matchPos[k1] <- elementMetadata(subSplit1)$posIndx[res1]

                  subSplit2 <- unlist(split(s1, strand(s1))[c("-")])
                  ## call follow
                  res2 <- follow(ranges(x1Split[[st]]), ranges(subSplit2))
                  indx2 <- !is.na(res2)
                  res2 <- res2[indx2]
                  k2 <- elementMetadata(x1Split[[st]])$posIndx[indx2]
                  matchPos[k2] <- elementMetadata(subSplit2)$posIndx[res2]
                  
                  mt <- k1[k1==k2]
                  for(p in mt) {
                       mn <- which.min( c(start(subSplit1[indx1])-end(ranges(x1Split[[st]])), 
                           start(ranges(x1Split[[st]])) - end(subSplit2[indx2])))
                       if(mn==1){
                           matchPos[p] <- elementMetadata(subSplit1)$posIndx[res1]
                       }
                  }
              }
          }
      }              
      matchPos
  })


setMethod("follow", c("GenomicRanges", "GenomicRanges"),
    function(x, subject, ...)
    {
        if(all(as.character(strand(x)) == "*") && all(as.character(strand(subject)) == "*"))
            strand(x) <- strand(subject) <- "+"
        
        elementMetadata(x) <- list("posIndx" = seq_len(length(x)))
        elementMetadata(subject) <- list("posIndx" = seq_len(length(subject)))
        xLst <- split(x, seqnames(x), drop = FALSE)
        xLst <- xLst[sapply(xLst, length) > 0]
        xSeqNames <- names(xLst)
        subjectLst <- split(subject, seqnames(subject), drop = FALSE)[xSeqNames]  
        matchPos <- rep.int(NA_integer_, length(x))

        for (sq in xSeqNames) {
            x1Split <- split(xLst[[sq]], strand(xLst[[sq]]))
            x1Split <- x1Split[lapply(x1Split, length) > 0]
            s1 <- subjectLst[[sq]]
            for (st in names(x1Split)) {
                if( st == "+" ){
                    subSplit <- unlist(split(s1, strand(s1))[c("+", "*")])
                    ## call follow
                    res <- follow(ranges(x1Split[[st]]), ranges(subSplit))
                    indx <- !is.na(res)
                    res <- res[indx]
                    matchPos[elementMetadata(x1Split[[st]])$posIndx[indx]] <-
                    elementMetadata(subSplit)$posIndx[res]

                }else if(st == "-"){
                    subSplit <- unlist(split(s1, strand(s1))[c("-", "*")])
                    ## call precede
                    res <- precede(ranges(x1Split[[st]]), ranges(subSplit))
                    indx <- !is.na(res)
                    res <- res[indx]
                    matchPos[elementMetadata(x1Split[[st]])$posIndx[indx]] <-
                    elementMetadata(subSplit)$posIndx[res]

                }else if(st == "*"){
                    subSplit1 <- unlist(split(s1, strand(s1))[c("+")])
                    ## call follow
                    res1 <- follow(ranges(x1Split[[st]]), ranges(subSplit1))
                    indx1 <- !is.na(res1)
                    res1 <- res1[indx1]
                    k1 <- elementMetadata(x1Split[[st]])$posIndx[indx1]
                    matchPos[k1] <- elementMetadata(subSplit1)$posIndx[res1]

                    subSplit2 <- unlist(split(s1, strand(s1))[c("-")])
                    ## call precede
                    res2 <- precede(ranges(x1Split[[st]]), ranges(subSplit2))
                    indx2 <- !is.na(res2)
                    res2 <- res2[indx2]
                    k2 <- elementMetadata(x1Split[[st]])$posIndx[indx2]
                    matchPos[k2] <- elementMetadata(subSplit2)$posIndx[res2]

                    mt <- k1[k1==k2]
                    for(p in mt) {
                       mn <- which.min( c(start(subSplit1[indx1])-end(ranges(x1Split[[st]])), 
                           start(ranges(x1Split[[st]])) - end(subSplit2[indx2])))
                       if(mn==1){
                           matchPos[p] <- elementMetadata(subSplit1)$posIndx[res1]
                       }
                    }
                }
            }
        }              
        matchPos
    })

