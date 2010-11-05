### =========================================================================
### Seqinfo objects
### -------------------------------------------------------------------------
###
### A Seqinfo object is a data.frame-like object that contains basic
### information about a set of genomic sequences. Currently only the
### length and circularity flag of each sequence is stored but more
### information might be added in the future.
###

setClass("Seqinfo",
    representation(
        seqnames="character",
        seqlengths="integer",
        is_circular="logical"
    )
)

### NOTE: Other possible (maybe better) names: GSeqinfo or GenomicSeqinfo
### for Genome/Genomic Sequence info.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setMethod("seqnames", "Seqinfo", function(x) x@seqnames)

setMethod("names", "Seqinfo", function(x) seqnames(x))

setMethod("length", "Seqinfo", function(x) length(seqnames(x)))

setMethod("seqlengths", "Seqinfo",
    function(x)
    {
        ans <- x@seqlengths
        names(ans) <- seqnames(x)
        ans
    }
)

### TODO: Put the isCircular(), isCircularWithKnownLength() and seqlengths()
### generics in the same file and same man page.
setGeneric("isCircular", function(x) standardGeneric("isCircular"))

setMethod("isCircular", "Seqinfo",
    function(x)
    {
        ans <- x@is_circular
        names(ans) <- seqnames(x)
        ans
    }
)

setGeneric("isCircularWithKnownLength",
    function(x) standardGeneric("isCircularWithKnownLength")
)

setMethod("isCircularWithKnownLength", "Seqinfo",
    function(x) ((isCircular(x) %in% TRUE) & !is.na(seqlengths(x)))
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.Seqinfo.seqnames <- function(x)
{
    x_seqnames <- seqnames(x)
    if (!is.character(x_seqnames)
     || !is.null(names(x_seqnames))
     || IRanges:::anyMissing(x_seqnames))
        return("'seqnames(x)' must be an unnamed character vector with no NAs")
    if (any(x_seqnames %in% "") || any(duplicated(x_seqnames)))
        return("'seqnames(x)' cannot contain zero-length or duplicated names")
    NULL
}

### Not really checking the slot itself but the value returned by the
### slot accessor.
.valid.Seqinfo.seqlengths <- function(x)
{
    x_seqlengths <- seqlengths(x)
    if (!is.integer(x_seqlengths)
     || length(x_seqlengths) != length(x)
     || !identical(names(x_seqlengths), seqnames(x)))
        return("'seqlengths(x)' must be an integer vector of the length of 'x' and with names 'seqnames(x)'")
    if (any(x_seqlengths < 0L, na.rm=TRUE))
        return("'seqlengths(x)' contains negative values")
    NULL
}

### Not really checking the slot itself but the value returned by the
### slot accessor.
.valid.Seqinfo.isCircular <- function(x)
{
    x_is_circular <- isCircular(x)
    if (!is.logical(x_is_circular)
     || length(x_is_circular) != length(x)
     || !identical(names(x_is_circular), seqnames(x)))
        return("'isCircular(x)' must be a logical vector of the length of 'x' and with names 'seqnames(x)'")
    NULL
}

.valid.Seqinfo <- function(x)
{
    c(.valid.Seqinfo.seqnames(x),
      .valid.Seqinfo.seqlengths(x),
      .valid.Seqinfo.isCircular(x))
}

setValidity2("Seqinfo", .valid.Seqinfo)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###
### Does only superficial checking of the arguments. The full validation is
### performed by new() thru the validity method.
###

### Make sure this always return an *unnamed* character vector.
.normargSeqnames <- function(seqnames)
{
    if (is.null(seqnames))
        return(character(0))
    if (!is.character(seqnames))
        stop("bad 'seqnames' value")
    unname(seqnames)
}

### Make sure this always return an *unnamed* integer vector.
.normargSeqlengths <- function(seqlengths, lx)
{
    if (identical(seqlengths, NA))
        return(rep.int(NA_integer_, lx))
    if (is.logical(seqlengths)) {
        if (all(is.na(seqlengths)))
            return(as.integer(seqlengths))
        stop("bad 'seqlengths' value")
    }
    if (!is.numeric(seqlengths))
        stop("bad 'seqlengths' value")
    if (!is.integer(seqlengths))
            return(as.integer(seqlengths))
    unname(seqlengths)
}

### Make sure this always return an *unnamed* logical vector.
.normargIsCircular <- function(isCircular, lx)
{
    if (identical(isCircular, NA))
        return(rep.int(NA, lx))
    if (!is.logical(isCircular))
        stop("bad 'isCircular' value")
    unname(isCircular)
}

Seqinfo <- function(seqnames=NULL, seqlengths=NA, isCircular=NA)
{
    seqnames <- .normargSeqnames(seqnames)
    seqlengths <- .normargSeqlengths(seqlengths, length(seqnames))
    is_circular <- .normargIsCircular(isCircular, length(seqnames))
    new("Seqinfo", seqnames=seqnames,
                   seqlengths=seqlengths,
                   is_circular=is_circular)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###

setReplaceMethod("seqnames", "Seqinfo",
    function(x, value)
    {
        value <- .normargSeqnames(value)
        if (length(value) != length(x))
            stop("'seqnames' replacement value ",
                 "must have 'length(x)' elements")
        x@seqnames <- value
        x
    }
)

setReplaceMethod("names", "Seqinfo",
    function(x, value) `seqnames<-`(x, value)
)

setReplaceMethod("seqlengths", "Seqinfo",
    function(x, value)
    {
        x@seqlengths <- .normargSeqlengths(value, length(x))
        x
    }
)

### TODO: Put the isCircular<-() and seqlengths<-() generics in the same
### file and same man page as the isCircular() and seqlengths() generics.
setGeneric("isCircular<-", function(x, value) standardGeneric("isCircular<-"))

setReplaceMethod("isCircular", "Seqinfo",
    function(x, value)
    {
        x@isCircular <- .normargIsCircular(value, length(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setMethod("as.data.frame", "Seqinfo",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        if (!is.null(row.names))
            warning("supplied 'row.names' value was ignored")
        if (!identical(optional, FALSE))
            warning("supplied 'optional' value was ignored")
        if (length(list(...)) != 0L)
            warning("extra arguments were ignored")
        data.frame(seqlengths=unname(seqlengths(x)),
                   isCircular=unname(isCircular(x)),
                   row.names=seqnames(x),
                   check.names=FALSE,
                   stringsAsFactors=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "show" method.
###

### cat(.showOutputAsCharacter(x), sep="\n") is equivalent to show(x).
.showOutputAsCharacter <- function(x)
{
    tmp <- tempfile()
    sink(file=tmp, type="output")
    show(x)
    sink(file=NULL)
    readLines(tmp)
}

.compactDataFrame <- function(x, max.nrow=19L)
{
    if (nrow(x) <= max.nrow)
        return(x)
    head <- head(x, n=max.nrow %/% 2L)
    tail <- tail(x, n=(max.nrow - 1L) %/% 2L)
    dotrow <- rep.int("...", ncol(x))
    names(dotrow) <- colnames(x)
    dotrow <- data.frame(as.list(dotrow),
                         row.names="...",
                         check.names=FALSE,
                         stringsAsFactors=FALSE)
    ## Won't handle properly the situation where one row in 'head' or 'tail'
    ## happens to be named "...".
    rbind(head, dotrow, tail)
}

### Should work properly on "narrow" data frames. Untested on data frames
### that are wider than the terminal.
showCompactDataFrame <- function(x, rownames.label="", left.margin="")
{
    compactdf <- .compactDataFrame(x)
    label_nchar <- nchar(rownames.label)
    if (label_nchar != 0L)
        row.names(compactdf) <- format(row.names(compactdf), width=label_nchar)
    showme <- .showOutputAsCharacter(compactdf)
    if (label_nchar != 0L)
        substr(showme[1L], 1L, label_nchar) <- rownames.label
    cat(paste(left.margin, showme, sep=""), sep="\n")
}

setMethod("show", "Seqinfo",
    function(object)
    {
        lo <- length(object)
        cat(class(object), " of length ", lo, "\n", sep="")
        if (lo == 0L)
            return(NULL)
        showCompactDataFrame(as.data.frame(object), "seqnames")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining.
###
### Why no "c" or "rbind" method? 'c(x, y)' would be expected to just append
### the rows in 'y' to the rows in 'x' resulting in an object of length
### 'length(x) + length(y)'. But that would not be compatible with the
### unicity of the "seqnames" key.
### So what we really need is a "merge" method that merge rows that have the
### same key in 'x' and 'y'.
###

.Seqinfo.mergelowhigh <- function(low, high)
{
    if (any(low != high, na.rm=TRUE))
        stop("incompatible Seqinfo objects")
    idx <- is.na(low) & !is.na(high)
    low[idx] <- high[idx]
    low
}

.Seqinfo.mergexy <- function(x, y)
{
    if (is.null(x)) {
        if (is.null(y))
            return(Seqinfo())
        if (is(y, "Seqinfo"))
            return(y)
        stop("all arguments must be Seqinfo objects (or NULLs)")
    }
    if (!is(x, "Seqinfo"))
        stop("all arguments must be Seqinfo objects (or NULLs)")
    if (is.null(y))
        return(x)
    if (!is(y, "Seqinfo"))
        stop("all arguments must be Seqinfo objects (or NULLs)")
    y2x_map <- match(seqnames(y), seqnames(x))
    ## Keep only rows from 'y' that are not already in 'x'.
    idx0 <- is.na(y2x_map)
    y0_seqnames <- seqnames(y)[idx0]
    y0_seqlengths <- seqlengths(y)[idx0]
    y0_is_circular <- isCircular(y)[idx0]
    ## Merge 'y' rows that already "exist" in 'x'.
    idx1 <- !idx0
    y2x_map1 <- y2x_map[idx1]
    x_seqlengths <- seqlengths(x)
    low <- x_seqlengths[y2x_map1]
    high <- seqlengths(y)[idx1]
    x_seqlengths[y2x_map1] <- .Seqinfo.mergelowhigh(low, high)
    x_is_circular <- isCircular(x)
    low <- x_is_circular[y2x_map1]
    high <- isCircular(y)[idx1]
    x_is_circular[y2x_map1] <- .Seqinfo.mergelowhigh(low, high)
    ## Make and return the result.
    Seqinfo(seqnames=c(seqnames(x), y0_seqnames),
            seqlengths=c(x_seqlengths, y0_seqlengths),
            isCircular=c(x_is_circular, y0_is_circular))
}

if (FALSE) {
.Seqinfo.merge <- function(...)
{
    args <- unname(list(...))
    ## Remove NULL elements...
    arg_is_null <- sapply(args, is.null)
    if (any(arg_is_null))
        args[arg_is_null] <- NULL  # ... by setting them to NULL!
    if (length(args) == 0L)
        return(Seqinfo())
    x <- args[[1L]]
    if (!all(sapply(args, is, class(x))))
        stop("all arguments in must be ", class(x), " objects (or NULLs)")

    if (length(args) == 1L)
        return(args[[1L]])
    ans
}

### These methods should not be called with named arguments: this tends to
### break dispatch!
setMethod("merge", c("Seqinfo", "missing"),
    function(x, y, ...) .Seqinfo.merge(x, ...)
)

setMethod("merge", c("Seqinfo", "Seqinfo"),
    function(x, y, ...) .Seqinfo.merge(x, y, ...)
)
}

