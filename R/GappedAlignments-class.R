### =========================================================================
### GappedAlignments objects
### -------------------------------------------------------------------------
###

setClass("GappedAlignments",
    contains="Sequence",
    representation(
        seqnames="Rle",               # 'factor' Rle
        start="integer",              # POS field in SAM
        cigar="character",            # extended CIGAR (see SAM format specs)
        strand="raw",
        seqinfo="Seqinfo"
        #mismatches="characterORNULL", # see MD optional field in SAM format specs
        #values="DataFrame"
    )
)

### Formal API:
###   length(x)   - single integer. Nb of alignments in 'x'.
###   rname(x)    - 'factor' Rle of the same length as 'x'.
###   seqnames(x) - same as 'rname(x)'.
###   rname(x) <- value - replacement form of 'rname(x)'.
###   seqnames(x) <- value - same as 'rname(x) <- value'.
###   cigar(x)    - character vector of the same length as 'x'.
###   strand(x)   - 'factor' Rle of the same length as 'x' (levels: +, -, *).
###   qwidth(x)   - integer vector of the same length as 'x'.
###   grglist(x)  - GRangesList object of the same length as 'x'.
###   granges(x), grg(x) - GRanges object of the same length as 'x'.
###   rglist(x)   - CompressedNormalIRangesList object of the same length as
###                 'x'.
###   ranges(x)   - IRanges object of the same length as 'x'.
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

setGeneric("granges", function(x) standardGeneric("granges"))

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

.GappedAlignmentsAsGRanges <- function(rname, start, width, strand, seqinfo)
{
    ranges <- IRanges(start=start, width=width)
    ans <- GRanges(seqnames=rname, ranges=ranges, strand=strand)
    seqinfo(ans) <- seqinfo
    ans
}

.GappedAlignmentsAsGRangesList <- function(rname, rglist, strand, seqinfo)
{
    nrg_per_alignment <- elementLengths(rglist)
    seqnames <- rep.int(rname, nrg_per_alignment)
    strand <- rep.int(strand, nrg_per_alignment)
    unlistData <- GRanges(seqnames=seqnames,
                          ranges=rglist@unlistData,
                          strand=strand)
    seqinfo(unlistData) <- seqinfo
    new("GRangesList",
        unlistData=unlistData,
        partitioning=rglist@partitioning,
        elementMetadata=new("DataFrame", nrows=length(rglist@partitioning)))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor-like methods.
###

setMethod("length", "GappedAlignments", function(x) length(x@cigar))

setMethod("rname", "GappedAlignments", function(x) x@seqnames)

### To make GappedAlignments/GRanges more interchangeable.
setMethod("seqnames", "GappedAlignments", function(x) rname(x))

setMethod("cigar", "GappedAlignments", function(x) x@cigar)

setMethod("width", "GappedAlignments", function(x) cigarToWidth(x@cigar))
setMethod("start", "GappedAlignments", function(x, ...) x@start)
setMethod("end", "GappedAlignments", function(x, ...) {x@start + width(x) - 1L})

setMethod("strand", "GappedAlignments",
    function(x)
        Rle(strand(IRanges:::compactBitvectorAsLogical(x@strand, length(x))))
)

setMethod("qwidth", "GappedAlignments", function(x) cigarToQWidth(x@cigar))

setMethod("grglist", "GappedAlignments",
    function(x)
        .GappedAlignmentsAsGRangesList(rname(x), rglist(x),
                                       strand(x), seqinfo(x))
)
setMethod("granges", "GappedAlignments",
    function(x)
        .GappedAlignmentsAsGRanges(rname(x), start(x), width(x),
                                   strand(x), seqinfo(x))
)

setMethod("grg", "GappedAlignments", function(x) granges(x))

setMethod("rglist", "GappedAlignments",
    function(x) cigarToIRangesListByAlignment(x@cigar, x@start)
)

setMethod("ranges", "GappedAlignments",
    function(x) IRanges(start=start(x), width=width(x))
)

setMethod("ngap", "GappedAlignments",
    function(x) {elementLengths(rglist(x)) - 1L}
)

setMethod("seqinfo", "GappedAlignments", function(x) x@seqinfo)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "rname<-" method.
###

.normargRNameReplaceValue <- function(x, value, ans.type=c("factor", "Rle"))
{
    ans.type <- match.arg(ans.type)
    if (!is.factor(value)
     && !is.character(value)
     && (!is(value, "Rle") || !is.character(runValue(value))
                              && !is.factor(runValue(value))))
        stop("'rname' value must be a character factor/vector, ",
             "or a 'character' Rle, or a 'factor' Rle")
    if (ans.type == "factor") {
        if (!is.factor(value))
            value <- as.factor(value)
    } else if (ans.type == "Rle") {
        ## We want to return a 'factor' Rle.
        if (!is(value, "Rle")) {
            if (!is.factor(value))
                value <- as.factor(value)
            value <- Rle(value)
        } else if (!is.factor(runValue(value))) {
            runValue(value) <- as.factor(runValue(value))
        }
    }
    if (length(value) != length(x))
        stop("'rname' value must be the same length as the object")
    value
}

### 'old_rname' and 'new_rname' must be 'factor' Rle.
.getRNameTranslationTable <- function(old_rname, new_rname)
{
    old <- runValue(old_rname)
    new <- runValue(new_rname)
    tmp <- unique(data.frame(old=old, new=new))
    if (!identical(runLength(old_rname), runLength(new_rname)) ||
        anyDuplicated(tmp$old) || anyDuplicated(tmp$new))
        stop("mapping between old an new 'rname' values is not one-to-one")
    if (isTRUE(all.equal(as.integer(tmp$old), as.integer(tmp$new)))) {
        tr_table <- levels(new)
        names(tr_table) <- levels(old)
    } else {
        tr_table <- tmp$new
        names(tr_table) <- tmp$old
    }
    tr_table
}

setReplaceMethod("rname", "GappedAlignments",
    function(x, value)
    {
        value <- .normargRNameReplaceValue(x, value, ans.type="Rle")
        tr_table <- .getRNameTranslationTable(rname(x), value)
        x@seqnames <- value
        seqnames(x@seqinfo) <- tr_table[seqlevels(x)]
        x
    }
)

### To make GappedAlignments/GRanges more interchangeable.
setReplaceMethod("seqnames", "GappedAlignments",
    function(x, value) `rname<-`(x, value)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other setters.
###

setReplaceMethod("seqinfo", "GappedAlignments",
    function(x, value)
    {
        if (!is(value, "Seqinfo"))
            stop("'value' must be a Seqinfo object")
        ## We only support this form of replacement for now.
        if (length(value) < length(seqinfo(x))
         || !identical(seqnames(value)[seq_len(length(seqinfo(x)))],
                       seqnames(seqinfo(x))))
            stop("the first elements in 'seqnames(value)' must be ",
                 "identical to 'seqnames(seqinfo(x))'")
        x@seqinfo <- value
        levels(x@seqnames) <- seqnames(value)
        validObject(x)
        x
    }
)

setReplaceMethod("seqlevels", "GappedAlignments",
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GappedAlignments.rname <- function(x)
{
    x_rname <- rname(x)
    if (!is(x_rname, "Rle") || !is.factor(runValue(x_rname))
     || !is.null(names(x_rname)) || any(is.na(x_rname)))
        return("'rname(x)' must be an unnamed 'factor' Rle with no NAs")
    if (length(x_rname) != length(cigar(x)))
        return("'rname(x)' and 'cigar(x)' must have the same length")
    NULL
}

.valid.GappedAlignments.start <- function(x)
{
    x_start <- start(x)
    if (!is.integer(x_start) || !is.null(names(x_start)) || IRanges:::anyMissing(x_start))
        return("'start(x)' must be an unnamed integer vector with no NAs")
    if (length(x_start) != length(cigar(x)))
        return("'start(x)' and 'cigar(x)' must have the same length")
    NULL
}

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

.valid.GappedAlignments.strand <- function(x)
{
    x_strand <- strand(x)
    if (!is(x_strand, "Rle") || !is.factor(runValue(x_strand))
     || !identical(levels(runValue(x_strand)), levels(strand()))
     || !is.null(names(x_strand)) || any(is.na(x_strand)))
        return("'strand(x)' must be an unnamed 'factor' Rle with no NAs (and with levels +, - and *)")
    if (length(x_strand) != length(cigar(x)))
        return("'strand(x)' and 'cigar(x)' must have the same length")
    NULL
}

.valid.GappedAlignments <- function(x)
{
    c(.valid.GappedAlignments.rname(x),
      .valid.GappedAlignments.start(x),
      .valid.GappedAlignments.cigar(x),
      .valid.GappedAlignments.strand(x),
      valid.GenomicRanges.seqinfo(x))
}

setValidity2("GappedAlignments", .valid.GappedAlignments,
             where=asNamespace("GenomicRanges"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

.asFactorRle <- function(x)
{
    if (is.character(x)) {
        x <- Rle(as.factor(x))
    } else if (is.factor(x)) {
        x <- Rle(x)
    } else if (is(x, "Rle") && is.character(runValue(x))) {
        runValue(x) <- as.factor(runValue(x))
    } else if (!is(x, "Rle") || !is.factor(runValue(x))) {
        stop("'x' must be a character vector, a factor, ",
             "a 'character' Rle, or a 'factor' Rle")
    }
    x
}

GappedAlignments <- function(rname=Rle(factor()), pos=integer(0),
                             cigar=character(0), strand=NULL, seqlengths=NULL)
{
    rname <- .asFactorRle(rname)
    if (any(is.na(rname)))
        stop("'rname' cannot have NAs")
    if (!is.integer(pos) || any(is.na(pos)))
        stop("'pos' must be an integer vector with no NAs")
    if (!is.character(cigar) || any(is.na(cigar)))
        stop("'cigar' must be a character vector with no NAs")
    if (is.null(strand)) {
        if (length(rname) != 0L)
            stop("'strand' must be specified when 'rname' is not empty")
        strand <- Rle(strand())
    } else if (is.factor(strand)) {
        strand <- Rle(strand)
    }
    strand <- IRanges:::logicalAsCompactBitvector(as.character(strand) == "-")
    if (is.null(seqlengths)) {
        seqlengths <- rep(NA_integer_, length(levels(rname)))
        names(seqlengths) <- levels(rname)
    } else if (!is.numeric(seqlengths)
            || is.null(names(seqlengths))
            || any(duplicated(names(seqlengths)))) {
        stop("'seqlengths' must be an integer vector with unique names")
    } else if (!setequal(names(seqlengths), levels(rname))) {
        stop("'names(seqlengths)' incompatible with 'levels(rname)'")
    } else if (!is.integer(seqlengths)) { 
        storage.mode(seqlengths) <- "integer"
    }
    seqinfo <- Seqinfo(seqnames=names(seqlengths), seqlengths=seqlengths)
    new("GappedAlignments", seqnames=rname, start=pos, cigar=cigar,
                            strand=strand, seqinfo=seqinfo)
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
### Coercion.
###

setAs("GappedAlignments", "GRangesList", function(from) grglist(from))
setAs("GappedAlignments", "GRanges", function(from) granges(from))
setAs("GappedAlignments", "RangesList", function(from) rglist(from))
setAs("GappedAlignments", "Ranges", function(from) ranges(from))


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
        if (is(i, "Rle"))
            i <- as.vector(i)
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
        x@seqnames <- x@seqnames[i]
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
        ans <- data.frame(rname=as.character(rname(x)),
                          strand=as.character(strand(x)),
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
              data.frame(rname=sketch(as.character(rname(object))),
                         strand=sketch(as.character(strand(object))),
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
        bounds <- try(callGeneric(granges(x), y), silent = TRUE)
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
