### =========================================================================
### GenomicRangesList objects
### -------------------------------------------------------------------------
###


setClass("GenomicRangesList",
    contains="RangesList",
    representation(
        "VIRTUAL",
        elementMetadata="DataFrame"
    ),
    prototype(
        elementType="GenomicRanges"
    )
)

.valid.GenomicRangesList.mcols <- function(x)
{
    msg <- NULL
    x_mcols <- x@elementMetadata
    if (nrow(x_mcols) != length(x))
        msg <- "'mcols(x)' has an incorrect number of rows"
    if (any(c("seqnames", "ranges", "strand", "start", "end", "width",
              "element") %in% colnames(x_mcols)))
        msg <-
          c(msg,
            paste("'mcols(x)' cannot have columns named \"seqnames\", ",
                  "\"ranges\", \"strand\", \"start\", \"end\", \"width\", ",
                  "or \"element\""))
    if (!is.null(rownames(x_mcols)))
        msg <- c(msg, "'mcols(x)' cannot have row names")
    msg
}

.valid.GenomicRangesList <- function(x)
{
    .valid.GenomicRangesList.mcols(x)
}

setValidity2("GenomicRangesList", .valid.GenomicRangesList)

### NOT exported but used in the GenomicAlignments package.
### FIXME: Seems to repeat a lot of code from IRanges:::show_IntegerRangesList!
show_GenomicRangesList <- function(x, with.header=TRUE)
{
    x_len <- length(x)
    if (with.header)
        cat(classNameForDisplay(x), " object of length ", x_len,
            ":\n", sep="")
    cumsumN <- end(PartitioningByEnd(x))
    N <- tail(cumsumN, 1)
    if (x_len == 0L) {
        cat("<0 elements>\n")
    } else if (x_len <= 3L || (x_len <= 5L && N <= 20L)) {
        ## Display full object.
        show(as.list(x))
    } else {
        ## Display truncated object.
        if (cumsumN[[3L]] <= 20L) {
            showK <- 3L
        } else if (cumsumN[[2L]] <= 20L) {
            showK <- 2L
        } else {
            showK <- 1L
        }
        show(as.list(x[seq_len(showK)]))
        diffK <- x_len - showK
        cat("...\n",
            "<", diffK, " more element",
            ifelse(diffK == 1L, "", "s"), ">\n",
            sep="")
    }
    #cat("-------\n")
    #cat("seqinfo: ", summary(seqinfo(x)), "\n", sep="")
}

setMethod("show", "GenomicRangesList",
    function(object) show_GenomicRangesList(object)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GenomicRanges_OR_GenomicRangesList
###

### TODO: Use this instead of GenomicRanges_OR_GRangesList in
### RangedSummarizedExperiment class definition to specify the class
### of its 'rowRanges' slot.
### FIXME: rtracklayer also defines GenomicRanges_OR_GenomicRangesList!
### Should use the class defined here instead.
setClassUnion("GenomicRanges_OR_GenomicRangesList",
    c("GenomicRanges", "GenomicRangesList")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SimpleGenomicRangesList and CompressedGenomicRangesList objects
###

setClass("SimpleGenomicRangesList",
    contains=c("GenomicRangesList", "SimpleRangesList"),
    representation("VIRTUAL")
)

setClass("CompressedGenomicRangesList",
    contains=c("GenomicRangesList", "CompressedRangesList"),
    representation("VIRTUAL")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GenomicRangesList getters
###

setMethod("seqnames", "GenomicRangesList",
    function(x)
    {
        if (!is(x, "CompressedList"))
            return(as(lapply(x, seqnames), "SimpleRleList"))
        unlisted_x <- unlist(x, use.names=FALSE)
        relist(seqnames(unlisted_x), x)
    }
)

setMethod("strand", "GenomicRangesList",
    function(x)
    {
        if (!is(x, "CompressedList"))
            return(as(lapply(x, strand), "SimpleRleList"))
        unlisted_x <- unlist(x, use.names=FALSE)
        relist(strand(unlisted_x), x)
    }
)

### NOT exported but used in the GenomicAlignments package.
### TODO: This should be defined as a method for RangesList objects (in
### the IRanges package). Then it would work on all RangesList derivatives,
### including IntegerRangesList, GenomicRangesList, and GAlignmentsList
### objects.
get_GenomicRangesList_mcols <-
    function(x, use.names=FALSE, level = c("between", "within"), ...)
{
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    level <- match.arg(level)
    if (level == "between") {
        ans <- x@elementMetadata
        if (use.names)
            rownames(ans) <- names(x)
        return(ans)
    }
    unlisted_x <- unlist(x, use.names=FALSE)
    unlisted_ans <- unlisted_x@elementMetadata
    if (use.names)
        rownames(unlisted_ans) <- names(unlisted_x)
    relist(unlisted_ans, x)
}
setMethod("elementMetadata", "GenomicRangesList", get_GenomicRangesList_mcols)

### seqinfo() is NOT supported on SimpleGenomicRangesList derivatives!
setMethod("seqinfo", "CompressedGenomicRangesList",
    function(x) seqinfo(x@unlistData)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GenomicRangesList setters
###
### Note that the setters are implemented on CompressedGenomicRangesList
### objects **only** at the moment. As a consequence, the setters do NOT
### work on SimpleGenomicRangesList objects!
###

### NOT exported but used in the GenomicAlignments package.
set_CompressedGenomicRangesList_seqnames <- function(x, value)
{
    if (!is(value, "AtomicList") ||
        !identical(elementNROWS(x), elementNROWS(value)))
        stop("replacement 'value' is not an AtomicList with the same ",
             "elementNROWS as 'x'")
    value <- unlist(value, use.names = FALSE)
    if (!is(value, "Rle"))
        value <- Rle(factor(value))
    else if (!is.factor(runValue(value)))
        runValue(value) <- factor(runValue(value))
    seqnames(x@unlistData) <- value
    x
}
setReplaceMethod("seqnames", "CompressedGenomicRangesList",
    set_CompressedGenomicRangesList_seqnames
)

### NOT exported but used in the GenomicAlignments package.
set_CompressedGenomicRangesList_strand <- function(x, value)
{
    if (!is(value, "AtomicList") ||
        !identical(elementNROWS(x), elementNROWS(value)))
        stop("replacement 'value' is not an AtomicList with the same ",
             "elementNROWS as 'x'")
    value <- unlist(value, use.names = FALSE)
    if (!is(value, "Rle"))
        value <- Rle(strand(value))
    else if (!is.factor(runValue(value)) ||
             !identical(levels(runValue(value)), levels(strand())))
        runValue(value) <- strand(runValue(value))
    strand(x@unlistData) <- value
    x
}
setReplaceMethod("strand", c("CompressedGenomicRangesList", "ANY"),
    set_CompressedGenomicRangesList_strand
)
setReplaceMethod("strand", c("CompressedGenomicRangesList", "character"),
    function(x, ..., value)
    {
        if (length(value) > 1L)
            stop("length(value) must be 1")
        strand(x@unlistData) <- value
        x
    }
)

### NOT exported but used in the GenomicAlignments package.
### TODO: This should be defined as a method for RangesList objects (in
### the IRanges package). Then it would work on all RangesList derivatives,
### including IntegerRangesList, GenomicRangesList, and GAlignmentsList
### objects.
set_CompressedGenomicRangesList_mcols <-
    function(x, level = c("between", "within"), ..., value)
{
    level <- match.arg(level)
    if (level == "between") {
        if (is.null(value))
            value <- S4Vectors:::make_zero_col_DataFrame(length(x))
        else if (!is(value, "DataFrame"))
            value <- DataFrame(value)
        if (!is.null(rownames(value)))
            rownames(value) <- NULL
        n <- length(x)
        k <- nrow(value)
        if (k != n) {
            if ((k == 0) || (k > n) || (n %% k != 0))
                stop(k, " rows in value to replace ", n, "rows")
            value <- value[rep(seq_len(k), length.out = n), , drop=FALSE]
        }
        x@elementMetadata <- value
    } else {
        if (is.null(value)) {
            value <- S4Vectors:::make_zero_col_DataFrame(length(x@unlistData))
        } else {
            if (!is(value, "SplitDataFrameList") ||
                !identical(elementNROWS(x), elementNROWS(value))) {
                stop("replacement 'value' is not a SplitDataFrameList with ",
                        "the same elementNROWS as 'x'")
            }
            value <- unlist(value, use.names = FALSE)
        }
        mcols(x@unlistData) <- value
    }
    x
}
setReplaceMethod("elementMetadata", "CompressedGenomicRangesList",
    set_CompressedGenomicRangesList_mcols
)

### NOT exported but used in the GenomicAlignments package.
set_CompressedGenomicRangesList_seqinfo <-
    function(x, new2old=NULL,
             pruning.mode=c("error", "coarse", "fine", "tidy"),
             value)
{
    pruning.mode <- match.arg(pruning.mode)
    if (!is(value, "Seqinfo"))
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    dangling_seqlevels <- GenomeInfoDb:::getDanglingSeqlevels(x,
                              new2old=new2old,
                              pruning.mode=pruning.mode,
                              seqlevels(value))
    if (length(dangling_seqlevels) != 0L) {
        ## Prune 'x'.
        non_dangling_range <- !(seqnames(x) %in% dangling_seqlevels)
        if (pruning.mode == "coarse") {
            x <- x[all(non_dangling_range)]
        } else {
            x <- x[non_dangling_range]  # "fine" pruning
            if (pruning.mode == "tidy") {
                ## Remove list elements that became empty because of "fine"
                ## pruning.
                x <- x[any(non_dangling_range) |
                       elementNROWS(non_dangling_range) == 0L]
            }
        }
    }
    seqinfo(x@unlistData, new2old=new2old) <- value
    x
}
setReplaceMethod("seqinfo", "CompressedGenomicRangesList",
    set_CompressedGenomicRangesList_seqinfo
)

### TODO: For the start, end, and width setters, we should probably have
### methods for SimpleRangesList and CompressedRangesList objects (in the
### IRanges package) instead of methods for IntegerRangesList and
### CompressedGenomicRangesList objects.

setReplaceMethod("start", "CompressedGenomicRangesList",
    function(x, ..., value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementNROWS(x), elementNROWS(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementNROWS as 'x'")
        value <- unlist(value, use.names=FALSE)
        start(x@unlistData, ...) <- value
        x
    }
)

setReplaceMethod("end", "CompressedGenomicRangesList",
    function(x, ..., value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementNROWS(x), elementNROWS(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementNROWS as 'x'")
        value <- unlist(value, use.names = FALSE)
        end(x@unlistData, ...) <- value
        x
    }
)

setReplaceMethod("width", "CompressedGenomicRangesList",
    function(x, ..., value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementNROWS(x), elementNROWS(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementNROWS as 'x'")
        value <- unlist(value, use.names = FALSE)
        width(x@unlistData, ...) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The score() accessor
###
### Kind of a silly accessor. And why is it defined at the GenomicRanges
### level and not at the Vector level? Or at least at the Ranges and
### RangesList levels so it works on IRanges and IRangesList objects too.
###

setMethod("score", "GenomicRangesList",
    function(x) mcols(x, use.names=FALSE)$score
)

setReplaceMethod("score", "GenomicRangesList",
    function(x, value)
    {
        mcols(x)$score <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods
###

### Do we still need the 2 methods below?

setMethod("updateObject", "CompressedGenomicRangesList",
    function(object, ..., verbose=FALSE)
    {
        ## unlistData slot.
        object@unlistData <- updateObject(object@unlistData,
                                          ..., verbose=verbose)
        ## Call method for CompressedList to update partitioning slot.
        callNextMethod()
    }
)

### Overwrite method for CompressedList objects just so we can fix on-the-fly
### any old GRanges instance (or other GenomicRanges derivative) stuck in
### the 'unlistData' slot.
setMethod("unlist", "CompressedGenomicRangesList",
    function(x, recursive=TRUE, use.names=TRUE)
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        unlisted_x <- updateObject(x@unlistData)
        if (use.names)
            unlisted_x <- S4Vectors:::set_unlisted_names(unlisted_x, x)
        unlisted_x
    }
)

