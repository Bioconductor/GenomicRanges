### =========================================================================
### Utility functions for subsetting and renaming seqlevels 
### -------------------------------------------------------------------------

keepSeqlevels <- function(x, value, ...)
{
    value <- unname(value)
    if (any(nomatch <- !value %in% seqlevels(x)))
        warning("invalid seqlevels '", 
                paste(value[nomatch], collapse=","), "' were ignored")
    seqlevels(x, force=TRUE) <- value[!nomatch]
    x
}

dropSeqlevels <- function(x, value, ...)
{
    value <- unname(value)
    if (any(nomatch <- !value %in% seqlevels(x)))
        warning("invalid seqlevels '", 
                paste(value[nomatch], collapse=","), "' were ignored")
    seqlevels(x, force=TRUE) <- seqlevels(x)[!seqlevels(x) %in% value]
    x
}

renameSeqlevels <- function(x, value, ...)
{
    nms <- names(value)
    ## unnamed
    if (is.null(nms)) {
        if (length(value) != length(seqlevels(x)))
            stop("unnamed 'value' must be the same length as seqlevels(x)")
        names(value) <- seqlevels(x)
    ## named
    } else {
        if (any(nomatch <- !nms %in% seqlevels(x)))
            warning("invalid seqlevels '", 
                    paste(nms[nomatch], collapse=","), "' were ignored")
        if (length(value) != length(seqlevels(x))) {
            level <- seqlevels(x)
            level[level %in% nms] <- value
            value <- level
        } 
    } 
    seqlevels(x) <- value  
    x 
}

## Currently applies to TranscriptDb only.
restoreSeqlevels <- function(x, ...)
{
    seqlevels0(x) 
    x
}

