### =========================================================================
### The seqinfo() (and releated) generic getters and setters
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqinfo() getter and setter.
###

setGeneric("seqinfo", function(x) standardGeneric("seqinfo"))

setGeneric("seqinfo<-", signature="x",
    function(x, old2new=NULL, value) standardGeneric("seqinfo<-")
)

### Compute the new seqnames associated with a seqinfo replacement.
### Assumes that 'seqnames(x)' is a 'factor' Rle (which is true if 'x' is a
### GRanges or GappedAlignments object, but not if it's a GRangesList object),
### and returns a 'factor' Rle of the same length (and same runLength vector).
makeNewSeqnames <- function(x, old2new=NULL, new_seqinfo)
{
    if (!is(new_seqinfo, "Seqinfo"))
        stop("supplied 'seqinfo' must be a Seqinfo object")
    x_seqnames <- seqnames(x)
    if (is.null(old2new)) {
        if (length(new_seqinfo) < length(seqinfo(x)) ||
            !identical(seqlevels(new_seqinfo)[seq_len(length(seqlevels(x)))],
                       seqlevels(x)))
            stop("when 'old2new' is NULL, the first sequence levels in the ",
                 "supplied 'seqinfo' must be identical to 'seqlevels(x)'")
        levels(x_seqnames) <- seqlevels(new_seqinfo)
        return(x_seqnames)
    }
    if (!is.integer(old2new) || length(old2new) != length(seqlevels(x)))
        stop("when 'old2new' is not NULL, it must be an integer ",
             "vector of the same length as 'seqlevels(x)'")
    min_old2new <- suppressWarnings(min(old2new, na.rm=TRUE))
    if (min_old2new != Inf) {
        if (min_old2new < 1L || max(old2new, na.rm=TRUE) > length(new_seqinfo))
            stop("non-NA values in 'old2new' must be >= 1 and <= N, ",
                 "where N is the nb of rows in the supplied 'seqinfo'")
    }
    if (any(duplicated(old2new) & !is.na(old2new)))
        stop("duplicates are not allowed among non-NA values in 'old2new'")
    dangling_seqlevels <- intersect(unique(x_seqnames),
                                    seqlevels(x)[is.na(old2new)])
    if (length(dangling_seqlevels) != 0L)
        stop("cannot drop levels currently in use: ",
             paste(dangling_seqlevels, collapse = ", "),
             ". Please consider subsetting 'x' first.")
    tmp <- runValue(x_seqnames)
    levels(tmp) <- seqlevels(new_seqinfo)[old2new]
    runValue(x_seqnames) <- factor(as.character(tmp),
                                   levels=seqlevels(new_seqinfo))
    return(x_seqnames)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqnames() getter and setter.
###

setGeneric("seqnames", function(x) standardGeneric("seqnames"))

setGeneric("seqnames<-", signature="x",
    function(x, value) standardGeneric("seqnames<-")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlevels() getter and setter.
###

setGeneric("seqlevels", function(x) standardGeneric("seqlevels"))

### Default "seqlevels" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("seqlevels", "ANY", function(x) seqlevels(seqinfo(x)))

setGeneric("seqlevels<-", signature="x",
    function(x, value) standardGeneric("seqlevels<-")
)

.is_partial_renaming <- function(old_seqlevels, new_seqlevels)
{
    if (length(new_seqlevels) != length(old_seqlevels))
        return(FALSE)
    is_renamed <- new_seqlevels != old_seqlevels
    tmp <- intersect(new_seqlevels[is_renamed], old_seqlevels[is_renamed])
    length(tmp) == 0L
}

### Returns -2 for "subsetting" mode, -1 for "renaming" mode, or an integer
### vector containing the mapping from the new to the old levels for "general"
### mode (i.e. a combination of renaming and/or subsetting). Note that this
### integer vector is guaranteed to contain no negative values.
getSeqlevelsReplacementMode <- function(seqlevels, old_seqlevels)
{
    ## Does NOT check for NA, duplicated or zero-length values since this will
    ## typically be done later by the Seqinfo() constructor.
    if (!is.character(seqlevels))
        stop("supplied 'seqlevels' must be a character vector")
    if (!is.null(names(seqlevels))) {
        nonempty_names <- names(seqlevels)[!(names(seqlevels) %in% c(NA, ""))]
        if (any(duplicated(nonempty_names)) ||
            length(setdiff(nonempty_names, old_seqlevels)) != 0L)
            stop("names of supplied 'seqlevels' contain duplicates ",
                 "or invalid sequence levels")
        return(match(names(seqlevels), old_seqlevels))
    }
    if (.is_partial_renaming(old_seqlevels, seqlevels))
        return(-1L)
    return(-2L)
}

### Returns the mapping from the old to the new sequence levels.
.getOld2NewSeqlevels <- function(old_seqlevels, new_seqlevels)
{
    mode <- getSeqlevelsReplacementMode(new_seqlevels, old_seqlevels)
    if (identical(mode, -2L))
        return(match(old_seqlevels, new_seqlevels))
    if (identical(mode, -1L))
        return(seq_len(length(new_seqlevels)))
    match(old_seqlevels, names(new_seqlevels))
}

### Default "seqlevels<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("seqlevels", "ANY",
    function(x, value)
    {
        x_seqinfo <- seqinfo(x)
        seqlevels(x_seqinfo) <- value
        old2new <- .getOld2NewSeqlevels(seqlevels(x), value)
        seqinfo(x, old2new=old2new) <- x_seqinfo
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqlengths() getter and setter.
###

setGeneric("seqlengths", function(x) standardGeneric("seqlengths"))

### Default "seqlengths" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("seqlengths", "ANY", function(x) seqlengths(seqinfo(x)))

setGeneric("seqlengths<-", signature="x",
    function(x, value) standardGeneric("seqlengths<-")
)

### Default "seqlengths<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("seqlengths", "ANY",
    function(x, value)
    {
        seqlengths(seqinfo(x)) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isCircular() getter and setter.
###

setGeneric("isCircular", function(x) standardGeneric("isCircular"))

### Default "isCircular" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("isCircular", "ANY", function(x) isCircular(seqinfo(x)))

setGeneric("isCircular<-", signature="x",
    function(x, value) standardGeneric("isCircular<-")
)

### Default "isCircular<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("isCircular", "ANY",
    function(x, value)
    {
        isCircular(seqinfo(x)) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### isCircularWithKnownLength() getter.
###
### TODO: Do we really need this?
###

setGeneric("isCircularWithKnownLength",
    function(x) standardGeneric("isCircularWithKnownLength")
)

### Default "isCircularWithKnownLength" method works on any object 'x' with
### a working "seqinfo" method.
setMethod("isCircularWithKnownLength", "ANY",
    function(x) isCircularWithKnownLength(seqinfo(x))
)

