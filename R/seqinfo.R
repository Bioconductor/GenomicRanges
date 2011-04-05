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

### Default "seqlevels<-" method works on any object 'x' with working
### "seqinfo" and "seqinfo<-" methods.
setReplaceMethod("seqlevels", "ANY",
    function(x, value)
    {
        ## More checkings of 'value' (e.g. no NAs, no zero-length or duplicated
        ## sequence names) is done later by the Seqinfo() constructor.
        if (!is.character(value))
            stop("supplied 'seqlevels' must be a character vector")
        if (!is.null(names(value))) {
            nonempty_names <- names(value)[!(names(value) %in% c(NA, ""))]
            if (any(duplicated(nonempty_names)) ||
                length(setdiff(nonempty_names, seqlevels(x))) != 0L)
                stop("names of supplied 'seqlevels' contain duplicates ",
                     "or invalid sequence levels")
        } else if (.is_partial_renaming(seqlevels(x), value)) {
            names(value) <- seqlevels(x)
        } else {
            ## Adding, dropping and/or reordering of the sequence levels
            ## *without* renaming.
            implicit_names <- value
            implicit_names[!(value %in% seqlevels(x))] <- ""
            names(value) <- implicit_names
        }
        old2new <- match(seqlevels(x), names(value))
        new2old <- rep.int(NA_integer_, length(value))
        new2old[old2new[!is.na(old2new)]] <-
            seq_len(length(seqlevels(x)))[!is.na(old2new)]
        new_seqlengths <- unname(seqlengths(x))[new2old]
        new_isCircular <- unname(isCircular(x))[new2old]
        new_seqinfo <- Seqinfo(value,
                               seqlengths=new_seqlengths,
                               isCircular=new_isCircular)
        seqinfo(x, old2new=old2new) <- new_seqinfo
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

