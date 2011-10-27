### =========================================================================
### Seqinfo objects
### -------------------------------------------------------------------------
###
### A Seqinfo object is a table-like object that contains basic information
### about a set of genomic sequences. The table has 1 row per sequence and
### 1 column per sequence attribute. Currently the only attributes are the
### length, circularity flag, and genome provenance (e.g. hg19) of the
### sequence, but more attributes might be added in the future as the need
### arises.
###

setClass("Seqinfo",
    representation(
        seqnames="character",
        seqlengths="integer",
        is_circular="logical",
        genome="character"
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

setMethod("seqlevels", "Seqinfo", function(x) seqnames(x))

setMethod("seqlengths", "Seqinfo",
    function(x)
    {
        ans <- x@seqlengths
        names(ans) <- seqnames(x)
        ans
    }
)

setMethod("isCircular", "Seqinfo",
    function(x)
    {
        ans <- x@is_circular
        names(ans) <- seqnames(x)
        ans
    }
)

setMethod("isCircularWithKnownLength", "Seqinfo",
    function(x) ((isCircular(x) %in% TRUE) & !is.na(seqlengths(x)))
)

setMethod("genome", "Seqinfo",
    function(x)
    {
        ans <- x@genome
        names(ans) <- seqnames(x)
        ans
    }
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

### Not really checking the slot itself but the value returned by the
### slot accessor.
.valid.Seqinfo.genome <- function(x)
{
    x_genome <- genome(x)
    if (!is.character(x_genome)
     || length(x_genome) != length(x)
     || !identical(names(x_genome), seqnames(x)))
        return("'genome(x)' must be a character vector of the length of 'x' and with names 'seqnames(x)'")
    NULL
}

.valid.Seqinfo <- function(x)
{
    c(.valid.Seqinfo.seqnames(x),
      .valid.Seqinfo.seqlengths(x),
      .valid.Seqinfo.isCircular(x),
      .valid.Seqinfo.genome(x))
}

setValidity2("Seqinfo", .valid.Seqinfo)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###
### The .normarg*() helper functions below do only partial checking of the
### arguments. The full validation is performed by new() thru the validity
### method.
###

### Make sure this always returns an *unnamed* character vector.
.normargSeqnames <- function(seqnames)
{
    if (is.null(seqnames))
        return(character(0))
    if (!is.character(seqnames))
        stop("bad supplied 'seqnames' vector")
    unname(seqnames)
}

### Make sure this always returns an *unnamed* integer vector.
.normargSeqlengths <- function(seqlengths, seqnames)
{
    if (identical(seqlengths, NA))
        return(rep.int(NA_integer_, length(seqnames)))
    if (!is.vector(seqlengths))
        stop("supplied 'seqlengths' must be a vector")
    if (length(seqlengths) != length(seqnames))
        stop("length of supplied 'seqlengths' must equal ",
             "the number of sequences")
    if (!is.null(names(seqlengths))
     && !identical(names(seqlengths), seqnames))
        stop("when the supplied 'seqlengths' vector is named, ",
             "the names must match the seqnames")
    if (is.logical(seqlengths)) {
        if (all(is.na(seqlengths)))
            return(as.integer(seqlengths))
        stop("bad supplied 'seqlengths' vector")
    }
    if (!is.numeric(seqlengths))
        stop("bad supplied 'seqlengths' vector")
    if (!is.integer(seqlengths))
        return(as.integer(seqlengths))
    unname(seqlengths)
}

### Make sure this always returns an *unnamed* logical vector.
.normargIsCircular <- function(isCircular, seqnames)
{
    if (identical(isCircular, NA))
        return(rep.int(NA, length(seqnames)))
    if (!is.vector(isCircular))
        stop("supplied 'isCircular' must be a vector")
    if (length(isCircular) != length(seqnames))
        stop("length of supplied 'isCircular' must equal ",
             "the number of sequences")
    if (!is.null(names(isCircular))
     && !identical(names(isCircular), seqnames))
        stop("when the supplied circularity flags are named, ",
             "the names must match the seqnames")
    if (!is.logical(isCircular))
        stop("bad supplied 'isCircular' vector")
    unname(isCircular)
}

### Make sure this always returns an *unnamed* character vector.
.normargGenome <- function(genome, seqnames)
{
    if (identical(genome, NA))
        return(rep.int(NA_character_, length(seqnames)))
    if (is.factor(genome))
        genome <- as.character(genome)
    else if (!is.vector(genome))
        stop("supplied 'genome' must be a vector")
    if (length(genome) != length(seqnames)) {
        if (length(genome) != 1L)
            stop("when length of supplied 'genome' vector is not 1, ",
                 "then it must equal the number of sequences")
        if (!is.null(names(genome)))
            stop("when length of supplied 'genome' vector is 1 ",
                 "and number of sequences is != 1, ",
                 "then 'genome' cannot be named")
        if (length(seqnames) == 0L)
            warning("supplied 'genome' vector has length 1 ",
                    "but number of sequences is 0")
        genome <- rep.int(genome, length(seqnames))
    } else if (!is.null(names(genome))
            && !identical(names(genome), seqnames))
        stop("when the supplied 'genome' vector is named, ",
             "the names must match the seqnames")
    if (is.logical(genome)) {
        if (all(is.na(genome)))
            return(as.character(genome))
        stop("bad supplied 'genome' vector")
    }
    if (!is.character(genome))
        stop("bad supplied 'genome' vector")
    unname(genome)
}

Seqinfo <- function(seqnames=NULL, seqlengths=NA, isCircular=NA, genome=NA)
{
    seqnames <- .normargSeqnames(seqnames)
    seqlengths <- .normargSeqlengths(seqlengths, seqnames)
    is_circular <- .normargIsCircular(isCircular, seqnames)
    genome <- .normargGenome(genome, seqnames)
    new("Seqinfo", seqnames=seqnames,
                   seqlengths=seqlengths,
                   is_circular=is_circular,
                   genome=genome)
}

setMethod("updateObject", "Seqinfo",
    function(object, ..., verbose=FALSE)
    {
        if (verbose)
            message("updateObject(object = 'Seqinfo')")
        if (!is(try(object@genome, silent=TRUE), "try-error"))
            return(genome)
        as(Seqinfo(seqnames=object@seqnames,
                   seqlengths=object@seqlengths,
                   isCircular=object@is_circular),
           class(object))
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###

### Support subsetting only by name.
setMethod("[", "Seqinfo",
    function(x, i, j, ..., drop=TRUE)
    {
        if (!missing(j) || length(list(...)) > 0L)
            stop("invalid subsetting")
        if (missing(i))
            return(x)
        if (!is.character(i))
            stop("a Seqinfo object can be subsetted only by name")
        if (!identical(drop, TRUE))
            warning("'drop' argument is ignored when subsetting ",
                    "a Seqinfo object")
        x_names <- names(x)
        i2names <- match(i, x_names)
        new_seqlengths <- unname(seqlengths(x))[i2names]
        new_isCircular <- unname(isCircular(x))[i2names]
        new_genome <- unname(genome(x))[i2names]
        Seqinfo(seqnames=i, seqlengths=new_seqlengths,
                isCircular=new_isCircular, genome=new_genome)

    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Setters.
###

setReplaceMethod("seqnames", "Seqinfo",
    function(x, value)
    {
        value <- .normargSeqnames(value)
        if (length(value) != length(x))
            stop("length of supplied 'seqnames' vector must equal ",
                 "the number of sequences")
        x@seqnames <- value
        x
    }
)

setReplaceMethod("names", "Seqinfo",
    function(x, value) `seqnames<-`(x, value)
)

setReplaceMethod("seqlevels", "Seqinfo",
    function(x, force=FALSE, value)
    {
        if (!identical(force, FALSE))
            warning("'force' is ignored in \"seqlevels<-\" method ",
                    "for Seqinfo objects")
        mode <- getSeqlevelsReplacementMode(value, seqlevels(x))
        if (identical(mode, -2L))
            return(x[value])
        if (identical(mode, -1L)) {
            seqnames(x) <- value
            return(x)
        }
        new_seqlengths <- unname(seqlengths(x))[mode]
        new_isCircular <- unname(isCircular(x))[mode]
        new_genome <- unname(genome(x))[mode]
        Seqinfo(seqnames=value, seqlengths=new_seqlengths,
                isCircular=new_isCircular, genome=new_genome)
    }
)

setReplaceMethod("seqlengths", "Seqinfo",
    function(x, value)
    {
        x@seqlengths <- .normargSeqlengths(value, seqnames(x))
        x
    }
)

setReplaceMethod("isCircular", "Seqinfo",
    function(x, value)
    {
        x@is_circular <- .normargIsCircular(value, seqnames(x))
        x
    }
)

setReplaceMethod("genome", "Seqinfo",
    function(x, value)
    {
        x@genome <- .normargGenome(value, seqnames(x))
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
                   genome=unname(genome(x)),
                   row.names=seqnames(x),
                   check.names=FALSE,
                   stringsAsFactors=FALSE)
    }
)

setAs("Seqinfo", "RangesList", function(from) {
  as(as(from, "GenomicRanges"), "RangesList")
})

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
### Why no "c" or "rbind" method for Seqinfo objects?
### 'c(x, y)' would be expected to just append the rows in 'y' to the rows in
### 'x' resulting in an object of length 'length(x) + length(y)'. But that
### would tend to break the constraint that the seqnames of a Seqinfo object
### must be unique.
### So what we really need is the ability to "merge" Seqinfo objects, that is,
### if a row in 'x' has the same seqname as a row in 'y', then the 2 rows must
### be merged in a single row before it's put in the result. If 2 rows cannot
### be merged because they contain incompatible information (e.g. different
### seqlengths or different circularity flags), then an error must be raised.
###

.Seqinfo.mergexy <- function(x, y)
{
    x_propernames <- setdiff(seqnames(x), seqnames(y))
    y_propernames <- setdiff(seqnames(y), seqnames(x))
    if (length(x_propernames) != 0L && length(y_propernames) != 0L) {
        msg <- c("Each of the 2 combined objects has sequence levels not in ",
                 "the other:\n  - in 'x': ",
                 paste(x_propernames, collapse=", "), "\n  - in 'y': ",
                 paste(y_propernames, collapse=", "), "\n",
                 "  Make sure to always combine/compare objects based on the ",
                 "same reference\n  genome (use suppressWarnings() to ",
                 "suppress this warning).")
        warning(msg)
    }
    ans_seqnames <- union(seqnames(x), seqnames(y))
    ans_seqlengths <- mergeNamedAtomicVectors(seqlengths(x), seqlengths(y),
                        what=c("sequence", "seqlengths"))
    ans_is_circular <- mergeNamedAtomicVectors(isCircular(x), isCircular(y),
                        what=c("sequence", "circularity flags"))
    ans_genome <- mergeNamedAtomicVectors(genome(x), genome(y),
                        what=c("sequence", "genomes"))
    Seqinfo(seqnames=ans_seqnames, seqlengths=ans_seqlengths,
            isCircular=ans_is_circular, genome=ans_genome)
}

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
    if (length(args) == 1L)
        return(x)
    args <- args[-1L]
    if (!all(sapply(args, is, class(x))))
        stop("all arguments in must be ", class(x), " objects (or NULLs)")
    for (y in args)
        x <- .Seqinfo.mergexy(x, y)
    x
}

### These methods should not be called with named arguments: this tends to
### break dispatch!
setMethod("merge", c("Seqinfo", "missing"),
    function(x, y, ...) .Seqinfo.merge(x, ...)
)

setMethod("merge", c("missing", "Seqinfo"),
    function(x, y, ...) .Seqinfo.merge(y, ...)
)

setMethod("merge", c("Seqinfo", "NULL"),
    function(x, y, ...) .Seqinfo.merge(x, ...)
)

setMethod("merge", c("NULL", "Seqinfo"),
    function(x, y, ...) .Seqinfo.merge(y, ...)
)

setMethod("merge", c("Seqinfo", "Seqinfo"),
    function(x, y, ...) .Seqinfo.merge(x, y, ...)
)

