### =========================================================================
### SeqInfo objects
### -------------------------------------------------------------------------
###
### A SeqInfo object is a data.frame-like object that contains basic
### information about a set of genomic sequences. Currently only the
### length and circularity flag of each sequence is stored but more
### information might be added in the future.
###

setClass("SeqInfo",
    representation(
        seqnames="character",
        seqlengths="integer",
        is_circular="logical"
    )
)

### NOTE: Other possible (maybe better) names: GSeqInfo or GenomicSeqInfo
### for Genome/Genomic Sequence Info.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters.
###

setMethod("seqnames", "SeqInfo", function(x) x@seqnames)

setMethod("names", "SeqInfo", function(x) seqnames(x))

setMethod("length", "SeqInfo", function(x) length(seqnames(x)))

setMethod("seqlengths", "SeqInfo",
    function(x)
    {
        ans <- x@seqlengths
        names(ans) <- seqnames(x)
        ans
    }
)

### TODO: Put the isCircular() and seqlengths() generics in the same file
### and same man page.
setGeneric("isCircular", function(x) standardGeneric("isCircular"))

setMethod("isCircular", "SeqInfo",
    function(x)
    {
        ans <- x@is_circular
        names(ans) <- seqnames(x)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.SeqInfo.seqnames <- function(x)
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
.valid.SeqInfo.seqlengths <- function(x)
{
    x_seqlengths <- seqlengths(x)
    if (!is.integer(x_seqlengths)
     || length(x_seqlengths) != length(x)
     || !identical(names(x_seqlengths), seqnames(x)))
        return("'seqlengths(x)' must be an integer vector of the ",
               "length of 'x' and with names 'seqnames(x)'")
    if (length(unique(is.na(x_seqlengths))) == 2L)
        return("'seqlengths(x)' cannot mix NA and non-NA values")
    if (any(x_seqlengths < 0L, na.rm=TRUE))
        return("'seqlengths(x)' contains negative values")
    NULL
}

### Not really checking the slot itself but the value returned by the
### slot accessor.
.valid.SeqInfo.isCircular <- function(x)
{
    x_is_circular <- isCircular(x)
    if (!is.logical(x_is_circular)
     || length(x_is_circular) != length(x)
     || !identical(names(x_is_circular), seqnames(x)))
        return("'isCircular(x)' must be a logical vector of the ",
               "length of 'x' and with names 'seqnames(x)'")
    if (length(unique(is.na(x_is_circular))) == 2L)
        return("'isCircular(x)' cannot mix NA and non-NA values")
    NULL
}

.valid.SeqInfo <- function(x)
{
    c(.valid.SeqInfo.seqnames(x),
      .valid.SeqInfo.seqlengths(x),
      .valid.SeqInfo.isCircular(x))
}

setValidity2("SeqInfo", .valid.SeqInfo)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###
### Does only superficial checking of the arguments. The full validation is
### performed by new() thru the validity method.
###

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

SeqInfo <- function(seqnames, seqlengths=NA, isCircular=NA)
{
    if (!is.character(seqnames))
        stop("'seqnames' must be a character vector")
    seqlengths <- .normargSeqlengths(seqlengths, length(seqnames))
    is_circular <- .normargIsCircular(isCircular, length(seqnames))
    new("SeqInfo", seqnames=seqnames,
                   seqlengths=seqlengths,
                   is_circular=is_circular)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###

setReplaceMethod("seqnames", "SeqInfo",
    function(x, value)
    {
        if (!is.vector(value) || !is.atomic(value))
            stop("bad 'seqnames' replacement value")
        value <- as.character(value)
        if (length(value) != length(x))
            stop("'seqnames' replacement value ",
                 "must have 'length(x)' elements")
        x@seqnames <- value
        x
    }
)

setReplaceMethod("names", "SeqInfo",
    function(x, value) `seqnames<-`(x, value)
)

setReplaceMethod("seqlengths", "SeqInfo",
    function(x, value)
    {
        x@seqlengths <- .normargSeqlengths(value, length(x))
        x
    }
)

### TODO: Put the isCircular<-() and seqlengths<-() generics in the same
### file and same man page as the isCircular() and seqlengths() generics.
setGeneric("isCircular<-", function(x, value) standardGeneric("isCircular<-"))

setReplaceMethod("isCircular", "SeqInfo",
    function(x, value)
    {
        x@isCircular <- .normargIsCircular(value, length(x))
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setMethod("as.data.frame", "SeqInfo",
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

setMethod("show", "SeqInfo",
    function(object)
    {
        lo <- length(object)
        cat(class(object), " of length ", lo, "\n", sep="")
        if (lo == 0L)
            return(NULL)
        showCompactDataFrame(as.data.frame(object), "seqnames")
    }
)

