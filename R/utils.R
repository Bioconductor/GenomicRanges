### =========================================================================
### Some low-level (non exported) utility functions.
### -------------------------------------------------------------------------


hasHead <- function(x, h) {
  identical(head(x, length(h)), h)
}

### TODO: Use this in GenomicFeatures::transcriptLocs2refLocs() and remove
### GenomicFeatures:::.normargExonStartsOrEnds().
### Used in GenomicAlignments package.
normargListOfIntegers <- function(arg, sep, argname) 
{
    if (is.list(arg))
        return(arg)
    if (is(arg, "IntegerList"))
        return(as.list(arg))
    if (is.character(arg))
        return(strsplitAsListOfIntegerVectors(arg, sep=sep))
    stop("'", argname, "' must be a list of integer vectors, ",
        "an IntegerList object,\n  or a character vector where ",
        "each element is a comma-separated list of\n  integers")
}

