### =========================================================================
### GappedAlignmentPairs objects
### -------------------------------------------------------------------------
###

### "first" and "last" GappedAlignments must have identical seqinfo.
setClass("GappedAlignmentPairs",
    contains="Vector",
    representation(
        NAMES="characterORNULL",      # R doesn't like @names !!
        first="GappedAlignments",     # of length N, no names, no elt metadata
        last="GappedAlignments",      # of length N, no names, no elt metadata
        isProperPair="logical",       # of length N
        elementMetadata="DataFrame"   # N rows
    )
)

### Formal API:
###   names(x)    - NULL or character vector.
###   length(x)   - single integer N. Nb of pairs in 'x'.
###   first(x)    - returns "first" slot.
###   last(x)     - returns "last" slot.
###   seqnames(x) - same as 'seqnames(first(x))' or 'seqnames(last(x))'.
###   strand(x)   - same as 'strand(first(x))' (opposite of 'strand(last(x))').
###   isProperPair(x) - returns "isProperPair" slot.
###   seqinfo(x)  - returns 'seqinfo(first(x))' (same as 'seqinfo(last(x))').
###   show(x)     - compact display in a data.frame-like fashion.
###   GappedAlignmentPairs(x) - constructor.
###   x[i]        - GappedAlignmentPairs object of the same class as 'x'
###                 (endomorphism).
###

setGeneric("first", function(x) standardGeneric("first"))
setGeneric("last", function(x) standardGeneric("last"))
setGeneric("isProperPair", function(x) standardGeneric("isProperPair"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setMethod("length", "GappedAlignmentPairs",
    function(x) length(x@first)
)

setMethod("names", "GappedAlignmentPairs",
    function(x) x@NAMES
)

setMethod("first", "GappedAlignmentPairs",
    function(x) setNames(x@first, names(x))
)

setMethod("last", "GappedAlignmentPairs",
    function(x) setNames(x@last, names(x))
)

setMethod("seqnames", "GappedAlignmentPairs",
    function(x) seqnames(first(x))
)

setMethod("strand", "GappedAlignmentPairs",
    function(x) strand(first(x))
)

setMethod("isProperPair", "GappedAlignmentPairs",
    function(x) x@isProperPair
)

setMethod("elementMetadata", "GappedAlignmentPairs",
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

setMethod("seqinfo", "GappedAlignmentPairs",
    function(x) seqinfo(x@first)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###

setReplaceMethod("names", "GappedAlignmentPairs",
    function(x, value)
    {
        if (!is.null(value))
            value <- as.character(value)
        x@NAMES <- value
        validObject(x)
        x
    }
)

setReplaceMethod("elementMetadata", "GappedAlignmentPairs",
    function(x, ..., value)
    {
        x@elementMetadata <- mk_elementMetadataReplacementValue(x, value)
        x
    }
)

setReplaceMethod("seqinfo", "GappedAlignmentPairs",
    function(x, new2old=NULL, force=FALSE, value)
    {
        if (!is(value, "Seqinfo"))
            stop("the supplied 'seqinfo' must be a Seqinfo object")
        first <- x@first
        last <- x@last
        dangling_seqlevels_in_first <- getDanglingSeqlevels(first,
                                           new2old=new2old, force=force,
                                           seqlevels(value))
        dangling_seqlevels_in_last <- getDanglingSeqlevels(last,
                                           new2old=new2old, force=force,
                                           seqlevels(value))
        dangling_seqlevels <- union(dangling_seqlevels_in_first,
                                    dangling_seqlevels_in_last)
        if (length(dangling_seqlevels) != 0L) {
            dropme_in_first <- seqnames(first) %in% dangling_seqlevels_in_first
            dropme_in_last <- seqnames(last) %in% dangling_seqlevels_in_last
            dropme <- dropme_in_first | dropme_in_last
            x <- x[!dropme]
        }
        seqinfo(x@first, new2old=new2old) <- value
        seqinfo(x@last, new2old=new2old) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GappedAlignmentPairs.names <- function(x)
{
    x_names <- names(x)
    if (is.null(x_names))
        return(NULL)
    if (!is.character(x_names) || !is.null(attributes(x_names))) {
        msg <- c("'names(x)' must be NULL or a character vector ",
                 "with no attributes")
        return(paste(msg, collapse=""))
    }
    if (length(x_names) != length(x))
        return("'names(x)' and 'x' must have the same length")
    NULL
}

.valid.GappedAlignmentPairs.first <- function(x)
{
    x_first <- x@first
    if (class(x_first) != "GappedAlignments")
        return("'x@first' must be a GappedAlignments instance")
    NULL
}

.valid.GappedAlignmentPairs.last <- function(x)
{
    x_last <- x@last
    if (class(x_last) != "GappedAlignments")
        return("'x@last' must be a GappedAlignments instance")
    x_first <- x@first
    if (length(x_last) != length(x_first))
        return("'x@last' and 'x@first' must have the same length")
    if (!identical(seqinfo(x_last), seqinfo(x_first)))
        return("'seqinfo(x@last)' and 'seqinfo(x@first)' must be identical")
    NULL
}

.valid.GappedAlignmentPairs.isProperPair <- function(x)
{
    x_isProperPair <- x@isProperPair
    if (!is.logical(x_isProperPair) || !is.null(attributes(x_isProperPair))) {
        msg <- c("'x@isProperPair' must be a logical vector ",
                 "with no attributes")
        return(paste(msg, collapse=""))
    }
    if (length(x_isProperPair) != length(x))
        return("'x@isProperPair' and 'x' must have the same length")
    if (IRanges:::anyMissing(x_isProperPair))
        return("'x@isProperPair' cannot contain NAs")
    NULL
}

.valid.GappedAlignmentPairs <- function(x)
{
    c(.valid.GappedAlignmentPairs.names(x),
      .valid.GappedAlignmentPairs.first(x),
      .valid.GappedAlignmentPairs.last(x),
      .valid.GappedAlignmentPairs.isProperPair(x))
}

setValidity2("GappedAlignmentPairs", .valid.GappedAlignmentPairs,
             where=asNamespace("GenomicRanges"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors.
###

GappedAlignmentPairs <- function(first, last, isProperPair, names=NULL)
{
    new2("GappedAlignmentPairs",
         NAMES=names,
         first=first, last=last,
         isProperPair=isProperPair,
         elementMetadata=new("DataFrame", nrows=length(first)),
         check=TRUE)
}

readGappedAlignmentPairs <- function(file, format="BAM", use.names=FALSE, ...)
{
    if (!isSingleString(file))
        stop("'file' must be a single string")
    if (!isSingleString(format))
        stop("'format' must be a single string")
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    if (format == "BAM") {
        suppressMessages(library("Rsamtools"))
        ans <- readBamGappedAlignmentPairs(file=file, use.names=use.names, ...)
        return(ans)
    }
    stop("only BAM format is supported at the moment")
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Supported 'i' types: numeric vector, logical vector, NULL and missing.
### TODO: Support subsetting by names.
setMethod("[", "GappedAlignmentPairs",
    function(x, i, j, ... , drop=TRUE)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        i <- IRanges:::normalizeSingleBracketSubscript(i, x)
        x@NAMES <- x@NAMES[i]
        x@first <- x@first[i]
        x@last <- x@last[i]
        x@isProperPair <- x@isProperPair[i]
        x@elementMetadata <- x@elementMetadata[i, , drop=FALSE]
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###

.makeNakedMatFromGappedAlignmentPairs <- function(x)
{
    lx <- length(x)
    nc <- ncol(elementMetadata(x))
    pair_cols <- cbind(seqnames=as.character(seqnames(x)),
                       strand=as.character(strand(x)))
    x_first <- first(x)
    first_cols <- cbind(ranges=IRanges:::showAsCell(ranges(x_first)))
    x_last <- last(x)
    last_cols <- cbind(ranges=IRanges:::showAsCell(ranges(x_last)))
    ans <- cbind(pair_cols,
                 `:`=rep.int(":", lx),
                 first_cols,
                 `:`=rep.int(":", lx),
                 last_cols)
    if (nc > 0L) {
        tmp <- do.call(data.frame, lapply(elementMetadata(x),
                                          IRanges:::showAsCell))
        ans <- cbind(ans, `|`=rep.int("|", lx), as.matrix(tmp))
    }
    ans
}

showGappedAlignmentPairs <- function(x, margin="",
                                        with.classinfo=FALSE,
                                        print.seqlengths=FALSE)
{
    lx <- length(x)
    nc <- ncol(elementMetadata(x))
    cat(class(x), " with ",
        lx, " alignment ", ifelse(lx == 1L, "pair", "pairs"),
        " and ",
        nc, " elementMetadata ", ifelse(nc == 1L, "col", "cols"),
        ":\n", sep="")
    out <- makePrettyMatrixForCompactPrinting(x,
               .makeNakedMatFromGappedAlignmentPairs)
    if (with.classinfo) {
        .PAIR_COL2CLASS <- c(
            seqnames="Rle",
            strand="Rle"
        )
        .HALVES_COL2CLASS <- c(
            ranges="IRanges"
        )
        .COL2CLASS <- c(.PAIR_COL2CLASS,
                        ":",
                        .HALVES_COL2CLASS,
                        ":",
                        .HALVES_COL2CLASS)
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

setMethod("show", "GappedAlignmentPairs",
    function(object)
        showGappedAlignmentPairs(object,
                                 with.classinfo=TRUE, print.seqlengths=TRUE)
)

