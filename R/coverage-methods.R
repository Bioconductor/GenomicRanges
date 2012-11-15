### =========================================================================
### coverage methods
### -------------------------------------------------------------------------
###

### Should work if 'arg' is an ordinary vector (i.e. atomic vector or list),
### or a List object.
.coverage.recycleAndSetNames <- function(arg, argname, seqlevels)
{
    k <- length(seqlevels)
    if (length(arg) < k)
        arg <- rep(arg, length.out = k)
    if (is.null(names(arg)))
        names(arg) <- seqlevels
    if (!all(seqlevels %in% names(arg)))
        stop("some seqnames missing from names(", argname, ")")
    arg
}

.coverage.normargShiftOrWeight <- function(arg, argname, x)
{
    if (is.list(arg) || is(arg, "List"))
        return(.coverage.recycleAndSetNames(arg, argname, seqlevels(x)))
    if (is.numeric(arg)) {
        if (length(arg) == 1L) {
            arg <- recycleNumericArg(arg, argname, length(seqlevels(x)))
            arg <- as.list(arg)
            names(arg) <- seqlevels(x)
        } else {
            arg <- recycleNumericArg(arg, argname, length(x))
            arg <- split(arg, as.factor(seqnames(x)))
        }
        return(arg)
    }
    stop("'", argname, "' must be a numeric vector, or a list, ",
         "or a List object")
}

.coverage.normargWidth <- function(width, x)
{
    if (is.null(width))
        width <- seqlengths(x)
    else if (!is.numeric(width))
        stop("'width' must be NULL or a numeric vector")
    .coverage.recycleAndSetNames(width, "width", seqlevels(x))
}

### 'circle.length' must be NA (if the underlying sequence is linear) or the
### length of the underlying circular sequence (integer vector of length 1
### with the name of the sequence).
### 'rg' must be an IRanges object.
.coverage.circle <- function(circle.length, rg, shift, width, weight)
{
    if (is.na(circle.length))
        return(coverage(rg, shift=shift, width=width, weight=weight))
    cvg <- fold(coverage(rg, weight=weight),
                circle.length, from=1L-shift)
    if (is.null(width))
        return(cvg)
    if (width > length(cvg))
        stop("invalid width (", width, ") ",
             "for circular sequence ", names(circle.length))
    cvg[seq_len(width)]
}

.coverage.seq <- function(seqlength, isCircular, rg, shift, width, weight)
{
    ## A shortcut for a common situation. Is it worth it?
    if (length(rg) == 0L && is.na(width))
        return(Rle(weight, 0L))
    if (isCircular %in% TRUE)
        circle.length <- seqlength
    else
        circle.length <- NA
    if (is.na(width))
        width <- NULL
    .coverage.circle(circle.length, rg, shift, width, weight)
}

setMethod("coverage", "GenomicRanges",
    function(x, shift=0L, width=NULL, weight=1L, ...)
    {
        if (any(start(x) < 1L))
            stop("'x' contains ranges starting before position 1. ",
                 "coverage() currently doesn't support this.")
        shift <- .coverage.normargShiftOrWeight(shift, "shift", x)
        width <- .coverage.normargWidth(width, x)
        if (isSingleString(weight)) {
            x_mcols <- mcols(x)
            if (!is(x_mcols, "DataTable")
             || sum(colnames(x_mcols) == weight) != 1L)
                stop("'mcols(x)' has 0 or more than 1 \"",
                     weight, "\" columns")
            weight <- x_mcols[[weight]]
        }
        weight <- .coverage.normargShiftOrWeight(weight, "weight", x)
        x_seqlevels <- seqlevels(x)
        x_seqlengths <- seqlengths(x)
        x_isCircular <- isCircular(x)
        x_ranges_by_seq <- split(unname(ranges(x)), seqnames(x))
        ii <- seq_len(length(x_seqlevels))
        names(ii) <- x_seqlevels
        ans_listData <- lapply(ii, function(i)
                               .coverage.seq(x_seqlengths[[i]],
                                             x_isCircular[[i]],
                                             x_ranges_by_seq[[i]],
                                             shift[[i]],
                                             width[[i]],
                                             weight[[i]]))
        IRanges:::newList("SimpleRleList", ans_listData)
    }
)

setMethod("coverage", "GRangesList",
    function(x, shift=0L, width=NULL, weight=1L, ...)
    {
        coverage(x@unlistData, shift=shift, width=width, weight=weight, ...)
    }
)

setMethod("coverage", "GappedAlignments",
    function(x, shift=0L, width=NULL, weight=1L, drop.D.ranges = FALSE, ...)
        coverage(grglist(x, drop.D.ranges = drop.D.ranges),
                 shift=shift, width=width, weight=weight, ...)
)

setMethod("coverage", "GappedAlignmentPairs",
    function(x, shift=0L, width=NULL, weight=1L, drop.D.ranges = FALSE, ...)
          coverage(grglist(x, drop.D.ranges = drop.D.ranges),
                   shift=shift, width=width, weight=weight, ...)
)

