### =========================================================================
### Utility functions for checking/normalizing user-supplied arguments
### -------------------------------------------------------------------------


### Return NULL or a Seqinfo object.
normarg_seqinfo1 <- function(seqinfo)
{
    if (is.null(seqinfo) || is(seqinfo, "Seqinfo"))
        return(seqinfo)
    if (is.character(seqinfo))
        return(Seqinfo(seqinfo))
    if (is.numeric(seqinfo)) {
        seqlevels <- names(seqinfo)
        if (is.null(seqlevels))
            stop(wmsg("when a numeric vector, 'seqinfo' must have names"))
        return(Seqinfo(seqlevels, seqlengths=seqinfo))
    }
    stop(wmsg("'seqinfo' must be NULL, or a Seqinfo object, ",
              "or a character vector of seqlevels, ",
              "or a named numeric vector of sequence lengths"))
}

### Return NULL or a Seqinfo object.
normarg_seqinfo2 <- function(seqinfo, seqlengths)
{
    seqinfo <- normarg_seqinfo1(seqinfo)
    if (is.null(seqlengths))
        return(seqinfo)
    ## Just a loose sanity check on 'seqlengths' before we call names() on
    ## it. The Seqinfo() constructor will take care of the full check and
    ## normalization by passing it thru GenomeInfoDb:::.normargSeqlengths().
    if (!is.vector(seqlengths))
        stop(wmsg("'seqlengths' must be NULL or a vector"))
    seqlengths_names <- names(seqlengths)
    if (is.null(seqlengths_names))
        stop(wmsg("'seqlengths' must have names"))
    seqinfo2 <- Seqinfo(seqlengths_names, seqlengths)
    if (is.null(seqinfo))
        return(seqinfo2)
    suppressWarnings(merge(seqinfo, seqinfo2))
}

### TODO: Use this in GenomicFeatures::transcriptLocs2refLocs() and
### remove GenomicFeatures:::.normargExonStartsOrEnds().
### Used in GenomicAlignments package.
normarg_list_of_integers <- function(arg, sep, argname)
{
    if (is.list(arg))
        return(arg)
    if (is(arg, "IntegerList"))
        return(as.list(arg))
    if (is.character(arg))
        return(toListOfIntegerVectors(arg, sep=sep))
    stop(wmsg("'", argname, "' must be a list of integer vectors, ",
              "an IntegerList object,\n  or a character vector where ",
              "each element is a comma-separated list of\n  integers"))
}

