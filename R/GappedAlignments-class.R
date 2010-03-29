### =========================================================================
### GappedAlignments objects
### -------------------------------------------------------------------------
###

setClass("GappedAlignments",
    contains="Sequence",
    representation(
        rname="factor",               # character factor
        start="integer",              # POS field in SAM
        cigar="character",            # extended CIGAR (see SAM format specs)
        strand="raw"
        #mismatches="characterORNULL", # see MD optional field in SAM format specs
        #values="DataFrame"
    )
)

### Formal API:
###   length(x)   - single integer. Nb of alignments in 'x'.
###   rname(x)    - character factor of the same length as 'x'.
###   rname(x) <- value - replacement form of 'rname(x)'.
###   strand(x)   - character factor of the same length as 'x' (levels: +, -,
###                 *).
###   cigar(x)    - character vector of the same length as 'x'.
###   qwidth(x)   - integer vector of the same length as 'x'.
###   grglist(x)  - GRangesList object of the same length as 'x'.
###   grg(x)      - GRanges object of the same length as 'x'.
###   rglist(x)   - CompressedNormalIRangesList object of the same length as
###                 'x'.
###   start(x), end(x), width(x) - integer vectors of the same length as 'x'.
###   ngap(x)     - integer vector of the same length as 'x'.
###   as.data.frame(x) - just a convenience used by show(x).
###   show(x)     - compact display in a data.frame-like fashion.
###   readGappedAlignments(x) - constructor.
###   x[i]        - GappedAlignments object of the same class as 'x'
###                 (endomorphism).
###
###   updateCigarAndStart(x, cigar=NULL, start=NULL) - GappedAlignments
###                 object of the same length and class as 'x' (endomorphism).
###                 For internal use only (NOT EXPORTED).
###
###   qnarrow(x, start=NA, end=NA, width=NA) - GappedAlignments object of the
###                 same length and class as 'x' (endomorphism).
###
###   narrow(x, start=NA, end=NA, width=NA) - GappedAlignments object of the
###                 same length and class as 'x' (endomorphism).
###
###   coverage(x) - named RleList object with one element (integer-Rle) per
###                 unique reference sequence.
###
###   findOverlaps(query, subject) - 'query' or 'subject' or both are
###                 GappedAlignments objects. Just a convenient wrapper for
###                 'findOverlaps(grglist(query), subject, ...)', etc...
###
###   countOverlaps(query, subject) - 'query' or 'subject' or both are
###                 GappedAlignments objects. Just a convenient wrapper for
###                 'countOverlaps(grglist(query), subject, ...)', etc...
###
###   subsetByOverlaps(query, subject) - 'query' or 'subject' or both are
###                 GappedAlignments objects.
###
### Concrete GappedAlignments implementations just need to implement:
###   length, rname, rname<-, strand, cigar, rglist, [ and updateCigarAndStart
### and the default methods defined in this file will work.

setGeneric("rname", function(x) standardGeneric("rname"))

setGeneric("rname<-", function(x, value) standardGeneric("rname<-"))

setGeneric("cigar", function(x) standardGeneric("cigar"))

setGeneric("qwidth", function(x) standardGeneric("qwidth"))

setGeneric("grglist", function(x) standardGeneric("grglist"))

setGeneric("grg", function(x) standardGeneric("grg"))

setGeneric("rglist", function(x) standardGeneric("rglist"))

setGeneric("updateCigarAndStart",
    function(x, cigar=NULL, start=NULL) standardGeneric("updateCigarAndStart")
)

setGeneric("qnarrow",
    function(x, start=NA, end=NA, width=NA) standardGeneric("qnarrow")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Two central helper functions.
###
### Note that their arguments are the different components of a
### GappedAlignments object instead of just the GappedAlignments object
### itself (arg 'x'). This allows them to be used in many different contexts
### e.g. when 'x' doesn't exist yet but is in the process of being constructed.
###

GappedAlignmentsAsGRanges <- function(rname, strand, start, width, check=TRUE)
{
    if (check) {
        if (is.factor(rname) && is.character(levels(rname)))
            rname <- as.character(rname)
        else if (!is.character(rname))
            stop("'rname' must be a character vector/factor")
        if (any(is.na(rname)))
            stop("'rname' cannot have NAs")
        if (!is.factor(strand) || !identical(levels(strand), levels(strand())))
            stop("'strand' must be a character factor")
        if (!is.integer(start))
            stop("'start' must be an integer vector")
        if (!is.integer(width))
            stop("'width' must be an integer vector")
    }
    seqnames <- Rle(rname)
    strand <- Rle(strand)
    ranges <- IRanges(start=start, width=width)
    GRanges(seqnames=seqnames, ranges=ranges, strand=strand)
}

GappedAlignmentsAsGRangesList <- function(rname, strand, rglist, check=TRUE)
{
    if (check) {
        if (is.factor(rname) && is.character(levels(rname)))
            rname <- as.character(rname)
        else if (!is.character(rname))
            stop("'rname' must be a character vector/factor")
        if (any(is.na(rname)))
            stop("'rname' cannot have NAs")
        if (!is.factor(strand) || !identical(levels(strand), levels(strand())))
            stop("'strand' must be a character factor")
        if (!is(rglist, "CompressedNormalIRangesList"))
            stop("'rglist' must be a CompressedNormalIRangesList object")
    }
    nrg_per_alignment <- elementLengths(rglist)
    seqnames <- Rle(rname, nrg_per_alignment)
    strand <- Rle(strand, nrg_per_alignment)
    unlistData <- GRanges(seqnames=seqnames,
                          ranges=rglist@unlistData,
                          strand=strand)
    new("GRangesList", unlistData=unlistData, partitioning=rglist@partitioning,
        elementMetadata = new("DataFrame", nrows = length(rglist@partitioning)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("length", "GappedAlignments", function(x) length(x@cigar))

setMethod("rname", "GappedAlignments", function(x) x@rname)

setMethod("cigar", "GappedAlignments", function(x) x@cigar)

setMethod("width", "GappedAlignments", function(x) cigarToWidth(x@cigar))
setMethod("start", "GappedAlignments", function(x, ...) x@start)
setMethod("end", "GappedAlignments", function(x, ...) {x@start + width(x) - 1L})

setMethod("strand", "GappedAlignments",
    function(x) strand(IRanges:::compactBitvectorAsLogical(x@strand, length(x)))
)

setMethod("qwidth", "GappedAlignments", function(x) cigarToQWidth(x@cigar))

setMethod("grglist", "GappedAlignments",
    function(x)
        GappedAlignmentsAsGRangesList(rname(x), strand(x), rglist(x))
)
setMethod("grg", "GappedAlignments",
    function(x)
        GappedAlignmentsAsGRanges(rname(x), strand(x), start(x), width(x))
)

setMethod("rglist", "GappedAlignments",
    function(x) cigarToIRangesListByAlignment(x@cigar, x@start)
)

setMethod("ngap", "GappedAlignments",
    function(x) {elementLengths(rglist(x)) - 1L}
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "rname<-" method.
###

normargRNameReplaceValue <- function(x, value, ans.type=c("factor", "Rle"))
{
    ans.type <- match.arg(ans.type)
    if (!is.factor(value)
     && !is.character(value)
     && (!is(value, "Rle") || !is.character(runValue(value))
                              && !is.factor(runValue(value))))
        stop("'rname' value must be a character factor/vector, ",
             "or a 'character' Rle, or a 'factor' Rle")
    if (ans.type == "factor" && !is.factor(value)) {
        value <- as.factor(value)
    } else if (ans.type == "Rle") {
        ## We want to return a 'character' Rle.
        if (!is(value, "Rle")) {
            if (!is.character(value))
                value <- as.character(value)
            value <- Rle(value)
        } else if (!is.character(runValue(value))) {
            runValue(value) <- as.character(runValue(value))
        }
    }
    if (length(value) != length(x))
        stop("'rname' value must be the same length as the object")

    ## Check the mapping between old and new values.
    old <- rname(x)  # guaranteed to be factor
    if (ans.type == "factor") {
        new <- value
    } else {
        old <- Rle(old)
        if (!identical(runLength(old), runLength(value))) {
            warning("mapping between old an new 'rname' values ",
                    "is not one-to-one")
            return(value)
        }
        old <- runValue(old)
        new <- runValue(value)
    }
    mapping <- unique(data.frame(old=old, new=new))
    if (any(duplicated(mapping$old)) || any(duplicated(mapping$new)))
        warning("mapping between old an new 'rname' values ",
                "is not one-to-one")
    value
}

setReplaceMethod("rname", "GappedAlignments",
    function(x, value)
    {
        x@rname <- normargRNameReplaceValue(x, value)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GappedAlignments.cigar <- function(x)
{
    x_cigar <- cigar(x)
    if (!is.character(x_cigar) || !is.null(names(x_cigar)) || any(is.na(x_cigar)))
        return("'cigar(x)' must be an unnamed character vector with no NAs")
    tmp <- validCigar(x_cigar)
    if (!is.null(tmp))
        return(paste("in 'cigar(x)':", tmp))
    NULL
}

.valid.GappedAlignments.rname <- function(x)
{
    x_rname <- rname(x)
    if (!is.factor(x_rname) || !is.character(levels(x_rname))
     || !is.null(names(x_rname)) || any(is.na(x_rname)))
        return("'rname(x)' must be an unnamed character factor with no NAs")
    if (length(x_rname) != length(cigar(x)))
        return("'rname(x)' and 'cigar(x)' must have the same length")
    NULL
}

.valid.GappedAlignments.strand <- function(x)
{
    x_strand <- strand(x)
    if (!is.factor(x_strand) || !identical(levels(x_strand), levels(strand()))
     || !is.null(names(x_strand)) || any(is.na(x_strand)))
        return("'strand(x)' must be an unnamed character factor with no NAs (and with levels +, - and *)")
    if (length(x_strand) != length(cigar(x)))
        return("'strand(x)' and 'cigar(x)' must have the same length")
    NULL
}

.valid.GappedAlignments.rglist <- function(x)
{
    x_rglist <- rglist(x)
    if (!is(x_rglist, "CompressedNormalIRangesList") || !is.null(names(x_rglist)))
        return("'rglist(x)' must be an unnamed CompressedNormalIRangesList object")
    if (length(x_rglist) != length(cigar(x)))
        return("'rglist(x)' and 'cigar(x)' must have the same length")
    if (any(elementLengths(x_rglist) == 0L))
        return("'rglist(x)' has elements with no ranges")
    x_start <- min(x_rglist)
    x_end <- max(x_rglist)
    if (!identical(x_end - x_start + 1L, width(x)))
        return("'rglist(x)' and 'width(x)' are out of sync")
    x_rglist2 <- cigarToIRangesListByAlignment(cigar(x), x_start)
    if (!identical(x_rglist2, x_rglist))
        return("'rglist(x)' and 'cigar(x)' are out of sync")
    NULL
}

.valid.GappedAlignments <- function(x)
{
    c(.valid.GappedAlignments.cigar(x),
      .valid.GappedAlignments.rname(x),
      .valid.GappedAlignments.strand(x),
      #.valid.GappedAlignments.grglist(x),
      #.valid.GappedAlignments.grg(x),
      .valid.GappedAlignments.rglist(x))
}

setValidity2("GappedAlignments", .valid.GappedAlignments,
             where=asNamespace("GenomicRanges"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

GappedAlignments <- function(rname=factor(), pos=integer(0),
                             cigar=character(0), strand=NULL)
{
    if (!is.factor(rname) || !is.character(levels(rname))) {
        if (!is.character(rname))
            stop("'rname' must be a character vector/factor")
        rname <- as.factor(rname)
    }
    if (any(is.na(rname)))
        stop("'rname' cannot have NAs")
    if (!is.integer(pos) || any(is.na(pos)))
        stop("'pos' must be an integer vector with no NAs")
    if (!is.character(cigar) || any(is.na(cigar)))
        stop("'cigar' must be a character vector with no NAs")
    if (is.null(strand)) {
        if (length(rname) != 0L)
            stop("'strand' must be specified when 'rname' is not empty")
        strand <- strand()
    } else if (!is.factor(strand)
           || !identical(levels(strand), levels(strand())))
        stop("invalid 'strand' argument")
    strand <- IRanges:::logicalAsCompactBitvector(strand == "-")
    new("GappedAlignments", rname=rname, start=pos, cigar=cigar, strand=strand)
}

readGappedAlignments <- function(file, format="BAM", ...)
{
    if (!isSingleString(file))
        stop("'file' must be a single string")
    if (!isSingleString(format))
        stop("'format' must be a single string")
    dotargs <- list(...)
    if (length(dotargs) != 0L && is.null(names(dotargs)))
        stop("extra arguments must be named")
    if (format == "BAM") {
        if ("index" %in% names(dotargs)) {
            index <- dotargs$index
            dotargs$index <- NULL
        } else {
            index <- file
        }
        if ("which" %in% names(dotargs)) {
            which <- dotargs$which
            dotargs$which <- NULL
        } else {
            which <- RangesList()
        }
        args <- c(list(file=file, index=index),
                  dotargs,
                  list(which=which))
        suppressMessages(library("Rsamtools"))
        ans <- do.call(readBamGappedAlignments, args)
        return(ans)
    }
    stop("only BAM format is supported for now")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
setMethod("[", "GappedAlignments",
    function(x, i, j, ... , drop=TRUE)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        if (!is.atomic(i))
            stop("invalid subscript type")
        lx <- length(x)
        if (length(i) == 0L) {
            i <- integer(0)
        } else if (is.numeric(i)) {
            if (min(i) < 0L)
                i <- seq_len(lx)[i]
            else if (!is.integer(i))
                i <- as.integer(i)
        } else if (is.logical(i)) {
            if (length(i) > lx)
                stop("subscript out of bounds")
            i <- seq_len(lx)[i]
        } else {
            stop("invalid subscript type")
        }
        x@rname <- x@rname[i]
        x@start <- x@start[i]
        x@cigar <- x@cigar[i]
        x@strand <- IRanges:::subsetCompactBitvector(x@strand, i)
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "as.data.frame" and "show" methods.
###

setMethod("as.data.frame", "GappedAlignments",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        if (!(is.null(row.names) || is.character(row.names)))
            stop("'row.names' must be NULL or a character vector")
        ans <- data.frame(rname=rname(x),
                          strand=strand(x),
                          cigar=cigar(x),
                          qwidth=qwidth(x),
                          start=start(x),
                          end=end(x),
                          width=width(x),
                          ngap=ngap(x),
                          row.names=row.names,
                          check.rows=TRUE,
                          check.names=FALSE,
                          stringsAsFactors=FALSE)
        return(ans)
    }
)

setMethod("show", "GappedAlignments",
    function(object)
    {
        lo <- length(object)
        cat(class(object), " of length ", lo, "\n", sep="")
        if (lo == 0L) {
            return(NULL)
        } else if (lo < 20L) {
            showme <-
              as.data.frame(object,
                            row.names=paste("[", seq_len(lo), "]", sep=""))
        } else {
            ## Use of as.vector() here is to prevent c() to do silly things
            ## when 'x' is a factor! (Try 'c(factor(LETTERS))', yes it's
            ## documented that c() will drop the attributes but still, this
            ## doesn't make sense).
            sketch <- function(x)
                          c(as.vector(x[1:9]),
                            "...",
                            as.vector(x[(length(x)-8L):length(x)]))
            showme <-
              data.frame(rname=sketch(rname(object)),
                         strand=sketch(strand(object)),
                         cigar=sketch(cigar(object)),
                         qwidth=sketch(qwidth(object)),
                         start=sketch(start(object)),
                         end=sketch(end(object)),
                         width=sketch(width(object)),
                         ngap=sketch(ngap(object)),
                         row.names=c(paste("[", 1:9, "]", sep=""), "...",
                                     paste("[", (lo-8L):lo, "]", sep="")),
                         check.rows=TRUE,
                         check.names=FALSE,
                         stringsAsFactors=FALSE)
        }
        show(showme)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "updateCigarAndStart" method.
###
### Performs atomic update of the cigar/start information.
### For internal use only (not exported).
###

setMethod("updateCigarAndStart", "GappedAlignments",
    function(x, cigar=NULL, start=NULL)
    {
        if (is.null(cigar))
            cigar <- cigar(x)
        else if (!is.character(cigar) || length(cigar) != length(x))
            stop("when not NULL, 'cigar' must be a character vector ",
                 "of the same length as 'x'")
        if (is.null(start))
            start <- start(x)
        else if (!is.integer(start) || length(start) != length(x))
            stop("when not NULL, 'start' must be an integer vector ",
                 "of the same length as 'x'")
        x@cigar <- cigar
        x@start <- start
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "qnarrow", "narrow", and "pintersect" methods.
###

setMethod("qnarrow", "GappedAlignments",
    function(x, start=NA, end=NA, width=NA)
    {
        ans_cigar <- cigarQNarrow(cigar(x),
                                  start=start, end=end, width=width)
        ans_start <- start(x) + attr(ans_cigar, "rshift")
        updateCigarAndStart(x, cigar=ans_cigar, start=ans_start)
    }
)

setMethod("narrow", "GappedAlignments",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
    {
        ans_cigar <- cigarNarrow(cigar(x),
                                 start=start, end=end, width=width)
        ans_start <- start(x) + attr(ans_cigar, "rshift")
        updateCigarAndStart(x, cigar=ans_cigar, start=ans_start)
    }
)

setMethod("pintersect", c("GappedAlignments", "GRanges"),
    function(x, y, ...)
    {
        bounds <- try(callGeneric(grg(x), y), silent = TRUE)
        if (inherits(bounds, "try-error"))
            stop("CIGAR is empty after intersection")
        start <- start(x) - (start(bounds) - 1L)
        start[which(start < 1L)] <- 1L
        end <- (end(bounds) - 1L) - end(x)
        end[which(end > -1L)] <- -1L
        narrow(x, start=start, end=end)
    }
)

setMethod("pintersect", c("GRanges", "GappedAlignments"),
    function(x, y, ...)
    {
        callGeneric(y, x)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "coverage" method.
###

setMethod("coverage", "GappedAlignments",
    function(x, shift=0L, width=NULL, weight=1L, ...)
    {
        shift <- as.list(shift)
        if (is.null(width))
            width <- list(width)
        else
            width <- as.list(width)
        weight <- as.list(weight)
        callGeneric(grglist(x), shift=shift, width=width, weight=weight, ...)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "findOverlaps" methods.
###

setMethod("findOverlaps", c("GappedAlignments", "ANY"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end"),
             select=c("all", "first"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        callGeneric(grglist(query), subject,
                    maxgap=maxgap, type=match.arg(type),
                    select=match.arg(select))
    }
)

setMethod("findOverlaps", c("ANY", "GappedAlignments"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end"),
             select=c("all", "first"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        callGeneric(query, grglist(subject),
                    maxgap=maxgap, type=match.arg(type),
                    select=match.arg(select))
    }
)

### Not strictly needed! Defining the above 2 methods covers that case but
### with the following note:
###   > findOverlaps(al1, al0)
###   Note: Method with signature "GappedAlignments#ANY" chosen for
###    function "findOverlaps", target signature
###    "GappedAlignments#GappedAlignments".
###    "ANY#GappedAlignments" would also be valid
setMethod("findOverlaps", c("GappedAlignments", "GappedAlignments"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end"),
             select=c("all", "first"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        callGeneric(grglist(query), grglist(subject),
                    maxgap=maxgap, type=match.arg(type),
                    select=match.arg(select))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "countOverlaps" methods.
###

setMethod("countOverlaps", c("GappedAlignments", "ANY"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        callGeneric(grglist(query), subject,
                    maxgap=maxgap, type=match.arg(type))
    }
)

setMethod("countOverlaps", c("ANY", "GappedAlignments"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        callGeneric(query, grglist(subject),
                    maxgap=maxgap, type=match.arg(type))
    }
)

setMethod("countOverlaps", c("GappedAlignments", "GappedAlignments"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        callGeneric(grglist(query), grglist(subject),
                    maxgap=maxgap, type=match.arg(type))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "subsetByOverlaps" methods.
###

setMethod("subsetByOverlaps", c("GappedAlignments", "ANY"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        type <- match.arg(type)
        query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                                  minoverlap = minoverlap, type = type,
                                  select = "first"))]
    }
)

setMethod("subsetByOverlaps", c("ANY", "GappedAlignments"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        type <- match.arg(type)
        query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                                  minoverlap = minoverlap, type = type,
                                  select = "first"))]
    }
)

setMethod("subsetByOverlaps", c("GappedAlignments", "GappedAlignments"),
    function(query, subject, maxgap=0L, minoverlap=1L,
             type=c("any", "start", "end"))
    {
        if (!identical(minoverlap, 1L))
            warning("'minoverlap' argument is ignored")
        type <- match.arg(type)
        query[!is.na(findOverlaps(query, subject, maxgap = maxgap,
                                  minoverlap = minoverlap, type = type,
                                  select = "first"))]
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "match" methods.
###

setMethod("match", c("GappedAlignments", "ANY"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL)
    {
        callGeneric(grglist(x), table,
                    nomatch=nomatch, incomparables=incomparables)
    }
)

setMethod("match", c("ANY", "GappedAlignments"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL)
    {
        callGeneric(x, grglist(table),
                    nomatch=nomatch, incomparables=incomparables)
    }
)

setMethod("match", c("GappedAlignments", "GappedAlignments"),
    function(x, table, nomatch = NA_integer_, incomparables = NULL)
    {
        callGeneric(grglist(x), grglist(table),
                    nomatch=nomatch, incomparables=incomparables)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "%in%" methods.
###

setMethod("%in%", c("GappedAlignments", "ANY"),
    function(x, table)
        callGeneric(grglist(x), table)
)

setMethod("%in%", c("ANY", "GappedAlignments"),
    function(x, table)
        callGeneric(x, grglist(table))
)

setMethod("%in%", c("GappedAlignments", "GappedAlignments"),
    function(x, table)
        callGeneric(grglist(x), grglist(table))
)
