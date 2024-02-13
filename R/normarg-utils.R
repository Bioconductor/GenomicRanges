### =========================================================================
### Utility functions for checking/normalizing user-supplied arguments
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### Used in the GenomicAlignments package!
### Return a factor-Rle with no NAs.
normarg_seqnames1 <- function(seqnames)
{
    if (is.null(seqnames))
        return(Rle(factor()))
    if (!is(seqnames, "Rle"))
        seqnames <- Rle(seqnames)
    run_vals <- runValue(seqnames)
    if (anyNA(run_vals))
        stop(wmsg("'seqnames' cannot contain NAs"))
    if (!is.factor(run_vals)) {
        if (!is.character(run_vals))
            run_vals <- as.character(run_vals)
        runValue(seqnames) <- factor(run_vals, levels=unique(run_vals))
    }
    seqnames
}

### Used in the GenomicAlignments package!
### 'seqnames' is assumed to be a factor-Rle with no NAs (which should
### be the case if it went thru normarg_seqnames1()).
### 'seqinfo' is assumed to be a Seqinfo object.
normarg_seqnames2 <- function(seqnames, seqinfo)
{
    ans_seqlevels <- seqlevels(seqinfo)
    run_vals <- runValue(seqnames)
    seqnames_levels <- levels(run_vals)
    is_used <- tabulate(run_vals, nbins=length(seqnames_levels)) != 0L
    seqnames_levels_in_use <- seqnames_levels[is_used]
    if (!all(seqnames_levels_in_use %in% ans_seqlevels))
        stop(wmsg("'seqnames' contains sequence names ",
                  "with no entries in 'seqinfo'"))
    if (!all(seqnames_levels %in% ans_seqlevels))
        warning(wmsg("levels in 'seqnames' with no entries ",
                     "in 'seqinfo' were dropped"))
    runValue(seqnames) <- factor(run_vals, levels=ans_seqlevels)
    seqnames
}

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

### Used in the GenomicAlignments package!
### Return NULL or a Seqinfo object.
normarg_seqinfo2 <- function(seqinfo, seqlengths)
{
    seqinfo <- normarg_seqinfo1(seqinfo)
    if (is.null(seqlengths))
        return(seqinfo)
    ## Just a loose sanity check on 'seqlengths' before we call names() on
    ## it. The Seqinfo() constructor will take care of the full check and
    ## normalization by passing it thru GenomeInfoDb:::.normarg_seqlengths().
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

### Used in the GenomicAlignments package!
### Return a factor-Rle with levels +|-|* and no NAs.
normarg_strand <- function(strand, ans_len)
{
    if (is.null(strand))
        return(Rle(strand("*"), ans_len))
    if (!is(strand, "Rle"))
        strand <- Rle(strand)
    run_vals <- runValue(strand)
    if (anyNA(run_vals)) {
        warning(wmsg("missing values in 'strand' converted to \"*\""))
        run_vals[is.na(run_vals)] <- "*"
    }
    if (!is.factor(run_vals) || !identical(levels(run_vals), levels(strand())))
        run_vals <- strand(run_vals)
    runValue(strand) <- run_vals
    strand
}

### Used in the GenomicAlignments package!
### TODO: Use this in GenomicFeatures::transcriptLocs2refLocs() and
### remove GenomicFeatures:::.normargExonStartsOrEnds().
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

