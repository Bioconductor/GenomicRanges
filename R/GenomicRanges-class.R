### =========================================================================
### The GenomicRanges interface
### -------------------------------------------------------------------------
###

### TODO: The 'constraint' slot could be moved to the Vector class (or to the
### Annotated class) so any Vector object could be constrained.
setClass("GenomicRanges",
    contains="Vector",
    representation(
        "VIRTUAL"#,
        #No more constraint slot for now...
        #constraint="ConstraintORNULL"
    )
)

setClassUnion("GenomicRangesORmissing", c("GenomicRanges", "missing"))

### The code in this file will work out-of-the-box on 'x' as long as
### seqnames(x), ranges(x), strand(x), seqlengths(x), seqinfo(),
### update(x) and clone(x) are defined.


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
    function(x, row.names=FALSE, ...)
    {
        if (!isTRUEorFALSE(row.names))
            stop("'row.names' must be TRUE or FALSE")
        ans <- x@elementMetadata
        if (row.names)
            rownames(ans) <- names(x)
        ans
    }
)

#setMethod("constraint", "GenomicRanges", function(x) x@constraint)


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

.valid.GenomicRanges.ranges <- function(x)
{
    if (class(ranges(x)) != "IRanges")
        return("'ranges(x)' must be an IRanges instance")
    NULL
}

.valid.GenomicRanges.strand <- function(x)
{
    if (!is.factor(runValue(strand(x))) ||
        !identical(levels(runValue(strand(x))), levels(strand())))
    {
        msg <- c("'strand' should be a 'factor' Rle with levels c(",
                 paste('"', levels(strand()), '"', sep="", collapse=", "),
                 ")")
        return(paste(msg, collapse=""))
    }
    if (IRanges:::anyMissing(runValue(strand(x))))
        return("'strand' contains missing values")
    NULL
}

## NOTE: This list is also included in the man page for GRanges objects.
## Keep the 2 lists in sync!
INVALID.GR.COLNAMES <- c("seqnames", "ranges", "strand",
                         "seqlevels", "seqlengths", "isCircular", "genome",
                         "start", "end", "width", "element")

.valid.GenomicRanges.elementMetadata <- function(x)
{    
    if (any(INVALID.GR.COLNAMES %in% colnames(elementMetadata(x)))) {
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
    if (any(x_seqlengths < 0L, na.rm=TRUE))
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
    if (any(minStarts[ncswkl] < 1L, na.rm=TRUE)
     || any(maxEnds[ncswkl] >
            seqlengths(x)[ncswkl], na.rm=TRUE))
        warning("'ranges' contains values outside of sequence bounds")
    NULL
}

.valid.GenomicRanges <- function(x)
{
    c(.valid.GenomicRanges.length(x),
      .valid.GenomicRanges.seqnames(x),
      .valid.GenomicRanges.ranges(x),
      .valid.GenomicRanges.strand(x),
      .valid.GenomicRanges.elementMetadata(x),
      valid.GenomicRanges.seqinfo(x))
      #checkConstraint(x, constraint(x)))
}

setValidity2("GenomicRanges", .valid.GenomicRanges)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("GenomicRanges", "RangedData",
    function(from)
    {
        rd <- RangedData(ranges(from), strand=strand(from),
                         elementMetadata(from), space=seqnames(from))
        elementMetadata(ranges(rd)) <- DataFrame(
                                         seqlengths=seqlengths(from),
                                         isCircular=isCircular(from),
                                         genome=genome(from))
        metadata(ranges(rd)) <- metadata(from)
        metadata(ranges(rd))$seqinfo <- seqinfo(from)
        rd
    }
)

setAs("GenomicRanges", "RangesList",
    function(from)
    {
        rl <- split(ranges(from), seqnames(from))
        emd <- split(DataFrame(strand=strand(from), elementMetadata(from)),
                     seqnames(from))
        rl <- mendoapply(function(ranges, metadata) {
          elementMetadata(ranges) <- metadata
          ranges
        }, rl, emd)
        elementMetadata(rl) <- DataFrame(seqlengths=seqlengths(from),
                                         isCircular=isCircular(from),
                                         genome=genome(from))
        metadata(rl) <- metadata(from)
        metadata(rl)$seqinfo <- seqinfo(from)
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
        data.frame(seqnames=as.factor(seqnames(x)),
                   start=start(x),
                   end=end(x),
                   width=width(x),
                   strand=as.factor(strand(x)),
                   as.data.frame(elementMetadata(x), ...),
                   row.names=row.names,
                   stringsAsFactors=FALSE)
    }
)

setAs("Seqinfo", "GenomicRanges",
    function(from)
    {
        if (any(is.na(seqlengths(from))))
            stop("cannot create a GenomicRanges from a Seqinfo ",
                 "with NA seqlengths")
        gr <- GRanges(seqnames(from),
                      IRanges(1L, width=seqlengths(from)),
                      seqlengths=seqlengths(from))
        seqinfo(gr) <- from
        names(gr) <- seqnames(gr)
        gr
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
        if (!identical(levels(value), seqlevels(x)))
            stop("levels of supplied 'seqnames' must be ",
                 "identical to 'seqlevels(x)'")
        n <- length(x)
        k <- length(value)
        if (k != n) {
            if ((k == 0L) || (k > n) || (n %% k != 0L))
                stop(k, " elements in value to replace ", n, " elements")
            value <- rep(value, length.out=n)
        }
        update(x, seqnames=value)
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
            if ((k == 0L) || (k > n) || (n %% k != 0L))
                stop(k, " elements in value to replace ", n, "elements")
            value <- rep(value, length.out=n)
        }
        update(x, ranges=value)
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
        strand <- rep(strand, length.out=n)
    }
    strand
}

setReplaceMethod("strand", "GenomicRanges",
    function(x, value) 
    {
        value <- normargGenomicRangesStrand(value, length(x))
        update(x, strand=value)
    }
)

setReplaceMethod("elementMetadata", "GenomicRanges",
    function(x, ..., value)
    {
        value <- mk_elementMetadataReplacementValue(x, value)
        update(x, elementMetadata=value)
    }
)

setReplaceMethod("seqinfo", "GenomicRanges",
    function(x, new2old=NULL, force=FALSE, value)
    {
        if (!is(value, "Seqinfo"))
            stop("the supplied 'seqinfo' must be a Seqinfo object")
        dangling_seqlevels <- getDanglingSeqlevels(x,
                                  new2old=new2old, force=force,
                                  seqlevels(value))
        if (length(dangling_seqlevels) != 0L)
            x <- x[!(seqnames(x) %in% dangling_seqlevels)]
        x <- update(x, seqnames=makeNewSeqnames(x, new2old, seqlevels(value)),
                       seqinfo=value)
        ## The ranges in 'x' need to be validated against
        ## the new sequence information (e.g. the sequence
        ## lengths might have changed).
        if (is.character(msg <- valid.GenomicRanges.seqinfo(x)))
          stop(msg)
        x
    }
)

setMethod("score", "GenomicRanges", function(x) values(x)$score)

#setReplaceMethod("constraint", "GenomicRanges",
#    function(x, value)
#    {
#        if (isSingleString(value))
#            value <- new(value)
#        if (!is(value, "ConstraintORNULL"))
#            stop("the supplied 'constraint' must be a ",
#                 "Constraint object, a single string, or NULL")
#        x@constraint <- value
#        validObject(x)
#        x
#    }
#)


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
    }
)

setGeneric("clone", function(x, ...) standardGeneric("clone"))

setMethod("clone", "ANY",
    function(x, ...)
    {
        if (nargs() > 1L)
            initialize(x, ...)
        else
            x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Ranges methods.
###

setMethod("start", "GenomicRanges", function(x, ...) start(ranges(x)))
setMethod("end", "GenomicRanges", function(x, ...) end(ranges(x)))
setMethod("width", "GenomicRanges", function(x) width(ranges(x)))

setReplaceMethod("start", "GenomicRanges",
    function(x, check=TRUE, value)
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
        start(ranges, check=check) <- starts
        update(x, ranges=ranges)
    }
)

setReplaceMethod("end", "GenomicRanges",
    function(x, check=TRUE, value)
    {
        if (!is.integer(value))
            value <- as.integer(value)
        ranges <- ranges(x)
        ends <- end(ranges)
        ends[] <- value
        seqlengths <- seqlengths(x)
        ## TODO: Revisit this to handle circularity.
        if (!IRanges:::anyMissing(seqlengths)) {
            seqlengths <- seqlengths[seqlevels(x)]
            maxEnds <- seqlengths[as.integer(seqnames(x))]
            trim <- which(ends > maxEnds)
            if (length(trim) > 0L) {
                warning("trimmed end values to be <= seqlengths")
                ends[trim] <- maxEnds[trim]
            }
        }
        end(ranges, check=check) <- ends
        update(x, ranges=ranges)
    }
)

setReplaceMethod("width", "GenomicRanges",
    function(x, check=TRUE, value)
    {
        if (!is.integer(value))
            value <- as.integer(value)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            end(x) <- start(x) + (value - 1L)
        } else {
            ranges <- ranges(x)
            width(ranges, check=check) <- value
            x <- update(x, ranges=ranges)
        }
        x
    }
)

setMethod("flank", "GenomicRanges",
    function(x, width, start=TRUE, both=FALSE, use.names=TRUE,
             ignore.strand=FALSE)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand)
            start <- rep.int(TRUE, length(x))
        else 
            start <- as.vector(start == (strand(x) != "-"))
        ranges <-
            flank(ranges(x), width=width, start=start, both=both,
                  use.names=use.names)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            start(x) <- start(ranges)
            end(x) <- end(ranges)
        } else {
            x <- clone(x, ranges=ranges)
        }
        x
    }
)

setMethod("resize", "GenomicRanges",
    function(x, width, fix="start", use.names=TRUE, ignore.strand=FALSE)
    {
        if (!missing(fix) &&
            (length(fix) > length(x) || length(x) %% length(fix) > 0L))
            stop("'x' is not a multiple of 'fix' length")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand) {
            fix <- Rle(rep.int(fix, length(x)))
        } else {
            revFix <- c(start="end", end="start", center="center")       
            fix <- ifelse(strand(x) == "-", revFix[fix], fix)
        }
        ranges <-
            resize(ranges(x), width=width, fix=fix, use.names=use.names)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            start(x) <- start(ranges)
            end(x) <- end(ranges)
        } else {
            x <- clone(x, ranges=ranges)
        }
        x
    }
)

## zooming (symmetrically scales the width)
setMethod("Ops", c("GenomicRanges", "numeric"),
    function(e1, e2)
    {
        if (IRanges:::anyMissing(e2))
            stop("NA not allowed as zoom factor")
        e2 <- recycleNumericArg(e2, "e2", length(e1))
        if (.Generic == "*") {
            e2 <- ifelse(e2 < 0, abs(1/e2), e2)
            resize(e1, width(e1) / e2, fix="center")
        } else {
            if (.Generic == "-") {
                e2 <- -e2
                .Generic <- "+"
            }
            if (.Generic == "+") {
                if (any(-e2*2 > width(e1)))
                    stop("adjustment would result in ranges ",
                         "with negative widths")
                resize(e1, width(e1) + e2*2, fix="center")
            }
        }
    }
)

setMethod("shift", "GenomicRanges",
    function(x, shift=0L, use.names=TRUE)
    {
        ranges <- shift(ranges(x), shift, use.names=use.names)
        if (!IRanges:::anyMissing(seqlengths(x))) {
            end(x, check=FALSE) <- end(ranges)
            start(x) <- pmin.int(start(ranges), end(x))
        } else {
            x <- clone(x, ranges=ranges)
        }
        x
    }
)

.interIntervalGenomicRanges <- function(x, FUN, ignore.strand=FALSE, ...)
{
    x <- clone(x)
    elementMetadata(x) <- NULL
    if (ignore.strand)
       f <- paste(seqnames(x), Rle(factor("*"), length(x)), sep="\r")
    else
       f <- paste(seqnames(x), strand(x), sep="\r")
    xIRangesList <- split(unname(ranges(x)), f)
    ansIRangesList <- FUN(xIRangesList, ...)
    k <- elementLengths(ansIRangesList)
    splitListNames <- strsplit(names(ansIRangesList), split="\r", fixed=TRUE)
    listNameMatrix <- matrix(as.character(unlist(splitListNames)), nrow=2L)
    ansSeqnames <- Rle(factor(listNameMatrix[1L, ], levels=seqlevels(x)), k)
    ansStrand <- Rle(strand(listNameMatrix[2L, ]), k)
    update(x, seqnames=ansSeqnames,
           ranges=unlist(ansIRangesList, use.names=FALSE),
           strand=ansStrand,
           elementMetadata=new("DataFrame", nrows=length(ansSeqnames)))
}

setMethod("disjoin", "GenomicRanges",
    function(x, ignore.strand=FALSE)
        .interIntervalGenomicRanges(x, disjoin, ignore.strand=ignore.strand)
)

applyOnRangesBySpace <- function(x, FUN, ..., ignore.strand = FALSE) {
  if (ignore.strand)
    f <- seqnames(x)
  else
    f <- paste(seqnames(x), strand(x), sep = "\r")
  xIRangesList <- split(unname(ranges(x)), f)
  ans <- FUN(xIRangesList, ...)
  if (is(ans, "List")) # values per range, otherwise assumed to be per-space
    ans <- unsplit(ans, f)
  ans
}

setMethod("isDisjoint", "GenomicRanges",
    function(x, ignore.strand=FALSE)
    {
        all(applyOnRangesBySpace(x, isDisjoint, ignore.strand = ignore.strand))
    }
)

setMethod("disjointBins", "GenomicRanges",
          function(x, ignore.strand = FALSE) {
            applyOnRangesBySpace(x, disjointBins,
                                 ignore.strand = ignore.strand)
          })

setMethod("gaps", "GenomicRanges",
    function(x, start=1L, end=seqlengths(x))
    {
        seqlevels <- seqlevels(x)
        if (!is.null(names(start)))
            start <- start[seqlevels]
        if (!is.null(names(end)))
            end <- end[seqlevels]
        start <- IRanges:::recycleVector(start, length(seqlevels))
        start <- rep(start, each=3L)
        end <- IRanges:::recycleVector(end, length(seqlevels))
        end <- rep(end, each=3L)
        .interIntervalGenomicRanges(x, gaps, start=start, end=end)
    }
)

setMethod("range", "GenomicRanges",
    function(x, ..., ignore.strand=FALSE, na.rm=FALSE)
        .interIntervalGenomicRanges(unname(c(x, ...)), range, 
                                    ignore.strand=ignore.strand)
)

setMethod("reduce", "GenomicRanges",
    function(x, drop.empty.ranges=FALSE, min.gapwidth=1L,
             with.inframe.attrib=FALSE, ignore.strand=FALSE)
    {
        if (!identical(with.inframe.attrib, FALSE))
            stop("'with.inframe.attrib' argument not supported ",
                 "when reducing a GenomicRanges object")
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (ignore.strand)
            strand(x) <- "*"
        .interIntervalGenomicRanges(x, reduce,
                                    drop.empty.ranges=drop.empty.ranges,
                                    min.gapwidth=min.gapwidth)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting and combining.
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
            if (!missing(j)) {
                elementMetadata <- elementMetadata(x, FALSE)[ , j, drop=FALSE]
                x <- clone(x, elementMetadata=elementMetadata)
            }
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
                if (length(whichEmpty) > 0L || !identical(nms, nms2))
                    names(ranges) <- nms2
            }
            x <- clone(x,
                       seqnames=seqnames(x)[i],
                       ranges=ranges,
                       strand=strand(x)[i],
                       elementMetadata=elementMetadata)
        }
        x
    }
)

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
                elementMetadata[ , ] <- elementMetadata(value, FALSE)
            else
                elementMetadata[ , j] <- elementMetadata(value, FALSE)
        } else {
            i <- iInfo[["idx"]]
            seqnames[i] <- seqnames(value)
            ranges[i] <- ranges(value)
            strand[i] <- strand(value)
            if (missing(j))
                elementMetadata[i, ] <- elementMetadata(value, FALSE)
            else
                elementMetadata[i, j] <- elementMetadata(value, FALSE)
        }
        update(x, seqnames=seqnames, ranges=ranges,
               strand=strand, elementMetadata=elementMetadata)
    }
)

setMethod("c", "GenomicRanges",
    function(x, ..., .ignoreElementMetadata=FALSE, recursive=FALSE)
    {
        if (recursive)
            stop("'recursive' mode not supported")
        args <- unname(list(x, ...))
        ans_seqinfo <- do.call(merge, lapply(args, seqinfo))
        ans_seqnames <- do.call(c, lapply(args, seqnames))
        ans_ranges <- do.call(c, lapply(args, ranges))
        ans_strand <- do.call(c, lapply(args, strand))
        
        if (.ignoreElementMetadata) {
          ans_elementMetadata <- new("DataFrame", nrows=length(ans_ranges))
        } else {
          ans_elementMetadata <- do.call(rbind,
                                         lapply(args, elementMetadata, FALSE))
        }
        
        ans_names <- names(ans_ranges)
        clone(x,
              seqnames=ans_seqnames,
              ranges=ans_ranges,
              strand=ans_strand,
              seqinfo=ans_seqinfo,
              elementMetadata=ans_elementMetadata)
    }
)

setMethod("seqselect", "GenomicRanges",
    function(x, start=NULL, end=NULL, width=NULL)
    {
        if (!is.null(end) || !is.null(width))
            start <- IRanges(start=start, end=end, width=width)
        irInfo <- IRanges:::.bracket.Index(start, length(x), names(x),
                                           asRanges=TRUE)
        if (!is.null(irInfo[["msg"]]))
            stop(irInfo[["msg"]])
        if (!irInfo[["useIdx"]])
            return(x)
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
        clone(x,
              seqnames=seqselect(seqnames(x), ir),
              ranges=ranges,
              strand=seqselect(strand(x), ir),
              elementMetadata=seqselect(elementMetadata(x, FALSE), ir))
    }
)

setReplaceMethod("seqselect", "GenomicRanges",
    function(x, start=NULL, end=NULL, width=NULL, value)
    {
        if (!is(value, "GenomicRanges"))
            stop("replacement value must be a GenomicRanges object")
        seqinfo(x) <- merge(seqinfo(x), seqinfo(value))
        if (is.null(end) && is.null(width)) {
            if (is.null(start))
                ir <- IRanges(start=1L, width=length(x))
            else if (is(start, "Ranges"))
                ir <- start
            else {
                if (is.logical(start) && length(start) != length(x))
                    start <- rep(start, length.out=length(x))
                ir <- as(start, "IRanges")
            }
        } else {
            ir <- IRanges(start=start, end=end, width=width, names=NULL)
        }
        ir <- reduce(ir)
        if (length(ir) == 0L) {
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
            update(x, seqnames=Rle(seqnames), ranges=ranges, 
                   strand=Rle(strand),
                   elementMetadata=elementMetadata)
        }
    }
)

setMethod("window", "GenomicRanges",
    function(x, start=NA, end=NA, width=NA,
             frequency=NULL, delta=NULL, ...)
    {
        update(x,
               seqnames=window(seqnames(x), start=start, end=end,
                               width=width, frequency=frequency,
                               delta=delta),
               ranges=window(ranges(x), start=start, end=end,
                             width=width, frequency=frequency,
                             delta=delta),
               strand=window(strand(x), start=start, end=end,
                             width=width, frequency=frequency,
                             delta=delta),
               elementMetadata=window(elementMetadata(x, FALSE),
                                      start=start, end=end,
                                      width=width, frequency=frequency,
                                      delta=delta))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###

.makeNakedMatFromGenomicRanges <- function(x)
{
    lx <- length(x)
    nc <- ncol(elementMetadata(x))
    ans <- cbind(seqnames=as.character(seqnames(x)),
                 ranges=IRanges:::showAsCell(ranges(x)),
                 strand=as.character(strand(x)))
    if (nc > 0L) {
        tmp <- do.call(data.frame, lapply(elementMetadata(x),
                                          IRanges:::showAsCell))
        ans <- cbind(ans, `|`=rep.int("|", lx), as.matrix(tmp))
    }
    ans
}

showGenomicRanges <- function(x, margin="",
                              with.classinfo=FALSE, print.seqlengths=FALSE)
{
    lx <- length(x)
    nc <- ncol(elementMetadata(x))
    cat(class(x), " with ",
        lx, " ", ifelse(lx == 1L, "range", "ranges"),
        " and ",
        nc, " elementMetadata ", ifelse(nc == 1L, "col", "cols"),
        ":\n", sep="")
    out <- makePrettyMatrixForCompactPrinting(x,
               .makeNakedMatFromGenomicRanges)
    if (with.classinfo) {
        .COL2CLASS <- c(
            seqnames="Rle",
            ranges="IRanges",
            strand="Rle"
        )
        classinfo <- makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
        ## A sanity check, but this should never happen!
        stopifnot(identical(colnames(classinfo), colnames(out)))
        out <- rbind(classinfo, out)
    }
    if (nrow(out) != 0L)
        rownames(out) <- paste(margin, rownames(out), sep="")
    print(out, quote=FALSE, right=TRUE)
    if (print.seqlengths) {
        cat(margin, "---\n", sep="")
        showSeqlengths(x, margin=margin)
    }
}

setMethod("show", "GenomicRanges",
    function(object)
        showGenomicRanges(object, margin="  ",
                          with.classinfo=TRUE, print.seqlengths=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other ranges utils.
###

setMethod("precede", c("GenomicRanges", "GenomicRanges"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        .findRanges(x, subject, ignore.strand, type="precede", ...) 
    }
)

setMethod("precede", c("GenomicRanges", "missing"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        callGeneric(x, subject=x, ignore.strand=ignore.strand, ...) 
    }
)

setMethod("follow", c("GenomicRanges", "GenomicRanges"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        .findRanges(x, subject, ignore.strand, type="follow", ...) 
    }
)

setMethod("follow", c("GenomicRanges", "missing"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        callGeneric(x, subject=x, ignore.strand=ignore.strand, ...) 
    }
)

setMethod("nearest", c("GenomicRanges", "GenomicRanges"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        .findRanges(x, subject, ignore.strand, type="nearest", ...) 
    }
)

setMethod("nearest", c("GenomicRanges", "missing"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        callGeneric(x, subject=x, ignore.strand=ignore.strand, ...)
    }
)

.findRanges <- function(x, subject, ignore.strand, type, ...)
{
    x <- .checkStrand(x, subject, ignore.strand)
    sinfo <- merge(seqinfo(x), seqinfo(subject))
    values(x) <- list(posIndx=seq_len(length(x)))
    values(subject) <- list(posIndx=seq_len(length(subject)))
    xLst <- split(x, seqnames(x), drop=FALSE)
    xLst <- xLst[elementLengths(xLst) > 0L]
    subjectLst <- split(subject, seqnames(subject), drop=FALSE)
    ans <- dist1 <- dist2 <- rep.int(NA_integer_, length(x))

    for (sq in names(xLst)) {
        if (!(sq %in% names(subjectLst)))
           break
        x1Split <- split(xLst[[sq]], strand(xLst[[sq]]))
        x1Split <- x1Split[elementLengths(x1Split) > 0L]
        s1 <- subjectLst[[sq]]
        ss1 <- split(s1, strand(s1))
        for (st in names(x1Split)) {
            xelt <- x1Split[[st]]
            if (st == "+") {
                subPos <- unlist(ss1[c("+", "*")])
                res <- .findRangesPos(xelt, subPos, type, ...)
                k <- values(xelt)[["posIndx"]][!is.na(res)]
                ans[k] <- values(subPos)[["posIndx"]][na.omit(res)]
            } else if (st == "-") {
                subNeg <- unlist(ss1[c("-", "*")])
                res <- .findRangesNeg(xelt, subNeg, type, ...)
                k <- values(xelt)[["posIndx"]][!is.na(res)]
                ans[k] <- values(subNeg)[["posIndx"]][na.omit(res)]
            } else if (st == "*") {
                if (type == "nearest") {
                    res <- nearest(ranges(xelt), ranges(s1))
                    k <- values(xelt)[["posIndx"]][!is.na(res)]
                    ans[k] <- values(s1)[["posIndx"]][na.omit(res)]
                } else { 
                    ## FIXME : include "*"?
                    sub1 <- unlist(ss1["+"])
                    res1 <- .findRangesPos(xelt, sub1, type, ...)
                    k1 <- values(xelt)[["posIndx"]][!is.na(res1)]
                    dist1[k1] <- .distance(xelt, sub1, res1, "+", type) 
                    ans[k1] <- values(sub1)[["posIndx"]][na.omit(res1)]

                    ## FIXME : include "*"?
                    sub2 <- unlist(ss1["-"])
                    res2 <- .findRangesNeg(xelt, sub2, type, ...)
                    k2 <- values(xelt)[["posIndx"]][!is.na(res2)]
                    dist2[k2] <- .distance(xelt, sub2, res2, "-", type) 
                    ans[k2] <- values(sub2)[["posIndx"]][na.omit(res2)]

                    ##  resolve queries matching both "+" and "-" 
                    if (!identical(integer(0), k1) && 
                        !identical(integer(0), k2)) {
                        tie <- .breakTie(sub1, res1, k1, dist1, sub2, res2, k2, 
                                         dist2, type)
                        ans[tie[["idx"]]] <- tie[["value"]]
                    }
                }
            }
        }
    }
    ans
}

.checkStrand <- function(x, subject, ignore.strand)
{
    if (!isTRUEorFALSE(ignore.strand))
        stop("'ignore.strand' must be TRUE or FALSE")
    if (ignore.strand)
        strand(x) <- strand(subject) <- "+"
    if (all(strand(x) == "*") && all(strand(subject) == "*"))
        strand(x) <- strand(subject) <- "+"
    x
}

.distance <- function(xelt, sub, res, strand, type)
{
    sidx <- na.omit(res)
    eidx <- !is.na(res)
    if (strand == "+")
        switch(type,
           precede=start(sub[sidx]) - end(xelt[eidx]),
           follow=start(xelt[eidx]) - end(sub[sidx]))
    else if (strand == "-") 
        switch(type,
           precede=start(xelt[eidx]) - end(sub[sidx]),
           follow=start(sub[sidx]) - end(xelt[eidx]))
}

.breakTie <- function(sub1, res1, k1, dist1, 
                      sub2, res2, k2, dist2, type)
{
    mt <- k1[k1 %in% k2]
    if (identical(integer(0), mt))
        return(DataFrame(idx=integer(), value=integer()))

    ## if tie, choose subject with minimum distance
    map1 <- data.frame(k1, spi=values(sub1)[["posIndx"]][na.omit(res1)])
    map2 <- data.frame(k2, spi=values(sub2)[["posIndx"]][na.omit(res2)])
    mindist <- mt[dist1[mt] < dist2[mt]]
    minvalue <- map1$spi[map1$k1 %in% mindist]

    ## if equidistant 
    ## precede -> choose subject with lowest index
    ## follow -> choose subject with highest index
    eqidist <- eqivalue <- integer()
    if (any(mt[dist1[mt] == dist2[mt]])) {
        eqidist <- mt[dist1[mt] == dist2[mt]]
        eq1 <- map1$spi[map1$k1 %in% eqidist]
        eq2 <- map2$spi[map2$k2 %in% eqidist]
        if (type == "precede")
            eqivalue <- pmin(eq1, eq2)
        if (type == "follow")
            eqivalue <- pmax(eq1, eq2)
    }
    DataFrame(idx=c(mindist, eqidist), value=c(minvalue, eqivalue)) 
}

.findRangesPos <- function(query, subject, type, ...)
{
    switch(type,
           precede=precede(ranges(query), ranges(subject), ...),
           follow=follow(ranges(query), ranges(subject), ...),
           nearest=nearest(ranges(query), ranges(subject), ...))
}

.findRangesNeg <- function(query, subject, type, ...)
{
    switch(type,
           precede=follow(ranges(query), ranges(subject), ...),
           follow=precede(ranges(query), ranges(subject), ...),
           nearest=nearest(ranges(query), ranges(subject), ...))
}

setMethod("narrow", "GenomicRanges",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        rng <- narrow(ranges(x), start=start, end=end, width=width,
                      use.names=TRUE)
        ranges(x) <- rng
        x
    }
)

.checkParms <- function(x, parm)
{
    if (!all(is.na(parm))) {
        if (!all(names(parm) %in% levels(seqnames(x))))
            stop("start should be a named numeric vector ",
                 "corresponding to seqnames")
    }
    temp <- structure(rep(NA_integer_, length(levels(seqnames(x)))), 
                      names=levels(seqnames(x)))
    temp[names(parm)] <- parm
    temp
}

.restrictRngs <- function(x, start, end, keep.all.ranges, use.names)
{
    tmp <- names(x)
    names(x) <- seq_len(length(x))
    rng <- ranges(x)
    res <- restrict(ranges(x), start, end, keep.all.ranges, use.names=TRUE)
    x <- x[as.numeric(names(res))]
    ranges(x) <- res
    if (!use.names)
        names(x) <- NULL
    else 
        names(x) <- tmp[as.numeric(names(res))]
    x
}

setMethod("restrict", "GenomicRanges",
    function(x, start=NA, end=NA, keep.all.ranges=FALSE, use.names=TRUE)
    {
        if (is.null(names(start)) && is.null(names(end)))
            return(.restrictRngs(x, start, end,keep.all.ranges, use.names))
        elementMetadata(x) <- cbind(elementMetadata(x),
                                    DataFrame(posIndx=seq_len(length(x))))
        splt <- split(x, seqnames(x))
        start <- .checkParms(x, start)
        end <- .checkParms(x, end) 
        res <- lapply(names(splt),
                      function(i) {
                          .restrictRngs(splt[[i]], start=start[i], end=end[i],
                                        keep.all.ranges, use.names)
                      })
        names(res) <- names(splt)
        ord <- unlist(GRangesList(res), use.names=FALSE)
        df <- data.frame(orig=elementMetadata(ord)$posIndx,
                         final=seq_len(length(ord)))
        indx <- with(df, order(orig, final))
        ord <- ord[indx, ]
        elementMetadata(ord) <- subset(elementMetadata(ord), select=-c(posIndx))
        ord
    }
)

setMethod("distance", c("GenomicRanges", "GenomicRanges"),
    function(x, y, ignore.strand=FALSE, ...)
    {
        if (!isTRUEorFALSE(ignore.strand))
            stop("'ignore.strand' must be TRUE or FALSE")
        if (length(x) != length(y))
            stop("'x' and 'y' must have the same length")
        d <- distance(ranges(x), ranges(y))
        mismtch <- as.character(seqnames(x)) != as.character(seqnames(y))
        if (any(mismtch))
            d[mismtch] <- NA
        if (!ignore.strand) {
            idx <- as.numeric(strand(x)) + as.numeric(strand(y))
            if (any(idx == 3))
                d[idx == 3] <- NA
        }
        d
    }
)

setMethod("distanceToNearest", c("GenomicRanges", "missing"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        callGeneric(x, subject=x, ignore.strand=ignore.strand, ...)
    }
)

setMethod("distanceToNearest", c("GenomicRanges", "GenomicRanges"),
    function(x, subject, ignore.strand=FALSE, ...)
    {
        x_nearest <- nearest(x, subject, ignore.strand=ignore.strand)
        if (identical(integer(0), x_nearest)) {
            DataFrame(queryHits=integer(), subjectHits=integer(),
                      distance=integer())
        } else {
            queryHits=seq(length(x))[!is.na(x_nearest)]
            subjectHits=na.omit(x_nearest)
            distance=distance(x[!is.na(x_nearest)], subject[na.omit(x_nearest)],
                              ignore.strand=ignore.strand)
            DataFrame(queryHits, subjectHits, distance)
        }
    }
)

