### =========================================================================
### "coverage" methods
### -------------------------------------------------------------------------
###

### Returns a list-like object.
.normarg_shift_or_weight <- function(arg, argname, x)
{
    if (is(arg, "list_OR_List")) {
        if (!identical(names(arg), seqlevels(x)))
            stop("when '", argname, "' is a list-like object, it must ",
                 "have 1 list element per seqlevel in 'x', and its names ",
                 "must be exactly 'seqlevels(x)'")
        return(arg)
    }
    arg <- IRanges:::replace_with_mcol_if_single_string(arg, x)
    if (!(is.numeric(arg) || is(arg, "Rle") && is.numeric(runValue(arg))))
        stop("'", argname, "' must be a numeric vector, a single string, ",
             "or a list-like object")
    arg <- S4Vectors:::V_recycle(arg, x, argname, "x")
    split(arg, seqnames(x))
}

setMethod("coverage", "GenomicRanges",
    function(x, shift=0L, width=NULL, weight=1L,
                method=c("auto", "sort", "hash", "naive"))
    {
        ## Normalize 'shift'.
        shift <- .normarg_shift_or_weight(shift, "shift", x)

        ## Just handle the default 'width' here. Non default will be checked
        ## in IRanges:::coverage_CompressedIRangesList().
        if (is.null(width))
            width <- seqlengths(x)

        ## Normalize 'weight'.
        weight <- .normarg_shift_or_weight(weight, "weight", x)

        x_ranges_list <- split(ranges(x), seqnames(x))
        circle.length <- seqlengths(x)
        circle.length[!(isCircular(x) %in% TRUE)] <- NA_integer_
        IRanges:::coverage_CompressedIRangesList(x_ranges_list,
                                        shift=shift,
                                        width=width,
                                        weight=weight,
                                        circle.length=circle.length,
                                        method=method,
                                        x_names.label="'seqlevels(x)'")
    }
)

### Overwrite above method with optimized method for StitchedGPos objects.
setMethod("coverage", "StitchedGPos",
    function(x, shift=0L, width=NULL, weight=1L,
                method=c("auto", "sort", "hash", "naive"))
    {
        CAN_ONLY_ETC <- c(" can only be a single number or a named ",
                          "list-like object when calling coverage() ",
                          "on a StitchedGPos object")
        if (!is(shift, "list_OR_List")) {
            if (!isSingleNumber(shift))
                stop(wmsg("'shift'", CAN_ONLY_ETC))
            shift <- Rle(shift)
        }
        if (!is(weight, "list_OR_List")) {
            if (!isSingleNumber(weight))
                stop(wmsg("'weight'", CAN_ONLY_ETC))
            weight <- Rle(weight)
        }
        x <- stitch_StitchedGPos(x)
        callGeneric()
    }
)

setMethod("coverage", "GRangesList",
    function(x, shift=0L, width=NULL, weight=1L,
                method=c("auto", "sort", "hash", "naive"))
    {
        coverage(x@unlistData, shift=shift, width=width, weight=weight,
                 method=method)
    }
)

