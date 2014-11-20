### =========================================================================
### Some low-level (non exported) utility functions.
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Used by "elementMetadata<-" methods.
###

### Try to turn 'value' into a DataFrame compatible with 'x'.
### Used in GenomicAlignments package.
normalizeMetadataColumnsReplacementValue <- function(value, x)
{
    if (is.null(value))
        return(new("DataFrame", nrows=length(x)))
    if (!is(value, "DataFrame"))
        value <- DataFrame(value)
    if (!is.null(rownames(value)))
        rownames(value) <- NULL
    n <- length(x)
    k <- nrow(value)
    if (k == n)
        return(value)
    if ((k == 0L) || (k > n) || (n %% k != 0L))
        stop(k, " rows in value to replace ", n, " rows")
    value[rep(seq_len(k), length.out=n), , drop=FALSE]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other stuff...
###

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

