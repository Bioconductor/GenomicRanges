### =========================================================================
### Some low-level (non exported) utility functions.
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Used by "show" methods.
###

makeClassinfoRowForCompactPrinting <- function(x, col2class)
{
    ans_names <- names(col2class)
    no_bracket <- ans_names == ""
    ans_names[no_bracket] <- col2class[no_bracket]
    left_brackets <- right_brackets <- character(length(col2class))
    left_brackets[!no_bracket] <- "<"
    right_brackets[!no_bracket] <- ">"
    ans <- paste0(left_brackets, col2class, right_brackets)
    names(ans) <- ans_names
    if (ncol(mcols(x)) > 0L) {
        tmp <- sapply(mcols(x),
                      function(xx) paste0("<", classNameForDisplay(xx), ">"))
        ans <- c(ans, `|`="|", tmp)
    }
    matrix(ans, nrow=1L, dimnames=list("", names(ans)))
}

compactPrintNamedAtomicVector <- function(x, margin="")
{
    x_len <- length(x)
    halfWidth <- (getOption("width") - nchar(margin)) %/% 2L
    first <- max(1L, halfWidth)
    showMatrix <-
      rbind(as.character(head(names(x), first)),
            as.character(head(x, first)))
    if (x_len > first) {
        last <- min(x_len - first, halfWidth)
        showMatrix <-
          cbind(showMatrix,
                rbind(as.character(tail(names(x), last)),
                      as.character(tail(x, last))))
    }
    showMatrix <- format(showMatrix, justify="right")
    cat(BiocGenerics:::labeledLine(margin, showMatrix[1L, ], count=FALSE,
                                           labelSep=""), sep="")
    cat(BiocGenerics:::labeledLine(margin, showMatrix[2L, ], count=FALSE,
                                           labelSep=""), sep="")
}

showSeqlengths <- function(object, margin="")
{
    cat(margin, "seqlengths:\n", sep="")
    compactPrintNamedAtomicVector(seqlengths(object), margin=margin)
}


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

