### =========================================================================
### The seqinfo() (and releated) generic getters and setters
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### seqinfo() getter and setter.
###

setGeneric("seqinfo", function(x) standardGeneric("seqinfo"))

setGeneric("seqinfo<-", signature="x",
    function(x, new2old=NULL, value) standardGeneric("seqinfo<-")
)

### Compute the new seqnames associated with a seqinfo replacement.
### Assumes that 'seqnames(x)' is a 'factor' Rle (which is true if 'x' is a
### GRanges or GappedAlignments object, but not if it's a GRangesList object),
### and returns a 'factor' Rle of the same length (and same runLength vector).
makeNewSeqnames <- function(x, new2old=NULL, new_seqinfo)
{
    if (!is(new_seqinfo, "Seqinfo"))
        stop("supplied 'seqinfo' must be a Seqinfo object")
    x_seqnames <- seqnames(x)
    M <- length(new_seqinfo)
    N <- length(seqlevels(x))
    if (is.null(new2old)) {
        if (M < N ||
            !identical(seqlevels(new_seqinfo)[seq_len(N)], seqlevels(x)))
            stop("when 'new2old' is NULL, the first sequence levels in the ",
                 "supplied 'seqinfo' must be indentical to 'seqlevels(x)'")
        levels(x_seqnames) <- seqlevels(new_seqinfo)
        return(x_seqnames)
    }
    if (!is.integer(new2old) || length(new2old) != M)
        stop("when 'new2old' is not NULL, it must be an integer ",
             "vector of the same length as the supplied 'seqinfo'")
    min_new2old <- suppressWarnings(min(new2old, na.rm=TRUE))
    if (min_new2old != Inf) {
        if (min_new2old < 1L || max(new2old, na.rm=TRUE) > N)
            stop("non-NA values in 'new2old' must be >= 1 and <= N, ",
                 "where N is the nb of sequence levels in 'x'")
    }
    if (any(duplicated(new2old) & !is.na(new2old)))
        stop("duplicates are not allowed among non-NA values in 'new2old'")
    old2new <- IRanges:::reverseIntegerInjection(new2old, N)
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

### Returns -2 for "subsetting" mode, -1 for "renaming" mode, or an integer
### vector containing the mapping from the new to the old levels for "general"
### mode (i.e. a combination of renaming and/or subsetting). Note that this
### integer vector is guaranteed to contain no negative values.
getSeqlevelsReplacementMode <- function(new_seqlevels, old_seqlevels)
{
    ## Does NOT check for NA, duplicated or zero-length values since this will
    ## typically be done later by the Seqinfo() constructor.
    if (!is.character(new_seqlevels))
        stop("supplied 'seqlevels' must be a character vector")
    nsl_names <- names(new_seqlevels)
    if (!is.null(nsl_names)) {
        nonempty_names <- nsl_names[!(nsl_names %in% c(NA, ""))]
        if (any(duplicated(nonempty_names)) ||
            length(setdiff(nonempty_names, old_seqlevels)) != 0L)
            stop("names of supplied 'seqlevels' contain duplicates ",
                 "or invalid sequence levels")
        return(match(nsl_names, old_seqlevels))
    }
    if (length(new_seqlevels) != length(old_seqlevels))
        return(-2L)
    is_renamed <- new_seqlevels != old_seqlevels
    tmp <- intersect(new_seqlevels[is_renamed], old_seqlevels[is_renamed])
    if (length(tmp) != 0L)
        return(-2L)
    return(-1L)
}

### Default "seqlevels<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("seqlevels", "ANY",
    function(x, value)
    {
        ## Make the new Seqinfo object.
        x_seqinfo <- seqinfo(x)
        seqlevels(x_seqinfo) <- value
        ## Map the new sequence levels to the old ones.
        new2old <- getSeqlevelsReplacementMode(value, seqlevels(x))
        if (identical(new2old, -2L)) {
            new2old <- match(value, seqlevels(x))
        } else if (identical(new2old, -1L)) {
            new2old <- seq_len(length(value))
        }
        ## Do the replacement.
        seqinfo(x, new2old=new2old) <- x_seqinfo
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### genome() getter and setter.
###

setGeneric("genome", function(x) standardGeneric("genome"))

### Default "genome" method works on any object 'x' with a working
### "seqinfo" method.
setMethod("genome", "ANY", function(x) genome(seqinfo(x)))

setGeneric("genome<-", signature="x",
    function(x, value) standardGeneric("genome<-")
)

### Default "genome<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("genome", "ANY",
    function(x, value)
    {
        genome(seqinfo(x)) <- value
        x
    }
)

