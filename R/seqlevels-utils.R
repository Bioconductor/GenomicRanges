### =========================================================================
### Utility functions for subsetting and renaming seqlevels 
### -------------------------------------------------------------------------

keepSeqlevels <- function(x, value, ...)
{
    value <- unname(value)
    if (any(nomatch <- !value %in% seqlevels(x)))
        warning("invalid seqlevels '", 
                paste(value[nomatch], collapse=","), "' were ignored")
    if (is(x, "BSgenome"))
        stop("seqlevels cannot be dropped from a BSgenome object")
    force <- is(x, "Seqinfo")
    seqlevels(x, force=!force) <- value[!nomatch]
    x
}

dropSeqlevels <- function(x, value, ...)
{
    value <- unname(value)
    if (any(nomatch <- !value %in% seqlevels(x)))
        warning("invalid seqlevels '", 
                paste(value[nomatch], collapse=","), "' were ignored")
    if (is(x, "BSgenome"))
        stop("seqlevels cannot be dropped from a BSgenome object")
    force <- is(x, "Seqinfo")
    seqlevels(x, force=!force) <- seqlevels(x)[!seqlevels(x) %in% value]
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

.checkStyleSpecies <- function(species, style, xlevels)
{
   require(GenomeInfoDb)
   #extract the standard chromsomes using species and style
   seqvec <- extractSeqlevels(species, style)
   seqvec[na.omit(match(xlevels,seqvec))]
}

keepStandardChromosomes <- function(x, species=species, style=style, ...)
{
    #check if user has entered same style as the given object
    originalstyle <- seqnameStyle(x)
    
    if(identical(originalstyle,style)){
        seq_vec <- .checkStyleSpecies(species=species,style=style,seqlevels(x))
        x <- keepSeqlevels(x,seq_vec)
        return(x)
    }else{
        stop("please check the seqnameStyle entered for object x.")
    }
}

