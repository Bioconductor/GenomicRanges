### =========================================================================
### coverage methods
### -------------------------------------------------------------------------
###

### TODO: Merge with Biostrings:::.V_recycle() and put in IRanges.
.recycle <- function(x, skeleton_len, x.label, skeleton.label)
{
    x_len <- length(x)
    if (x_len == skeleton_len)
        return(x)
    if (x_len < skeleton_len) {
        if (x_len == 0L)
            stop("cannot recycle zero-length '", x.label, "' ",
                 "to the length of '", skeleton.label, "'")
    } else {
        if (x_len >= 2L)
            stop("'", x.label, "' is longer than '", skeleton.label, "'")
    }
    if (skeleton_len %% x_len != 0L)
        warning("'", x.label, "' length is not a divisor ",
                "of '", skeleton.label, "' length")
    rep(x, length.out=skeleton_len)
}

.normarg_shift_or_weight <- function(arg, arg.label, x)
{
    if (is.list(arg) || is(arg, "List")) {
        if (!identical(names(arg), seqlevels(x)))
            stop("when '", arg.label, "' is a list-like object, it must ",
                 "have 1 list element per seqlevel in 'x', and its names ",
                 "must be exactly 'seqlevels(x)'")
        return(arg)
    }
    if (isSingleString(arg)) {
        x_mcols <- mcols(x)
        if (!is(x_mcols, "DataTable")
         || sum(colnames(x_mcols) == arg) != 1L)
            stop("'mcols(x)' has 0 or more than 1 \"",
                 arg, "\" columns")
        arg <- x_mcols[ , arg]
    }
    if (!is.numeric(arg))
        stop("'", arg.label, "' must be a numeric vector, a single string, ", 
             "or a list-like object")
    split(.recycle(arg, length(x), arg.label, "x"), seqnames(x))
}

setMethod("coverage", "GenomicRanges",
    function(x, shift=0L, width=NULL, weight=1L,
                method=c("auto", "sort", "hash"))
    {
        ## Normalize 'shift'.
        shift <- .normarg_shift_or_weight(shift, "shift", x)

        ## Normalize 'width'.
        if (is.null(width)) {
            width <- seqlengths(x)
        } else if (!(is.numeric(width)
                  && identical(names(width), seqlevels(x)))) {
            stop("'width' must be NULL or an integer vector with ",
                 "the length and names of 'seqlengths(x)'")
        }

        ## Normalize 'weight'.
        weight <- .normarg_shift_or_weight(weight, "weight", x)

        x_ranges_list <- split(ranges(x), seqnames(x))
        circle.length <- seqlengths(x)
        circle.length[!(isCircular(x) %in% TRUE)] <- NA_integer_
        IRanges:::.CompressedIRangesList.coverage(x_ranges_list,
                                        shift=shift,
                                        width=width,
                                        weight=weight,
                                        circle.length=circle.length,
                                        method=method)
    }
)

setMethod("coverage", "GRangesList",
    function(x, shift=0L, width=NULL, weight=1L,
                method=c("auto", "sort", "hash"))
    {
        coverage(x@unlistData, shift=shift, width=width, weight=weight,
                 method=method)
    }
)

setMethod("coverage", "GAlignments",
    function(x, shift=0L, width=NULL, weight=1L,
                method=c("auto", "sort", "hash"), drop.D.ranges=FALSE)
        coverage(grglist(x, drop.D.ranges=drop.D.ranges),
                 shift=shift, width=width, weight=weight, method=method)
)

setMethod("coverage", "GAlignmentPairs",
    function(x, shift=0L, width=NULL, weight=1L,
                method=c("auto", "sort", "hash"), drop.D.ranges=FALSE)
        coverage(grglist(x, drop.D.ranges=drop.D.ranges),
                 shift=shift, width=width, weight=weight, method=method)
)

