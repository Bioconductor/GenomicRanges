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
                 paste("\"", INVALID.GR.COLNAMES, "\"", sep="", collapse=", "),
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
        tmp <- do.call(data.frame, c(lapply(elementMetadata(x),
                                            IRanges:::showAsCell),
                                     list(check.names=FALSE)))
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

