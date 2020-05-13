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

OLD_GRANGESLIST_INSTANCE_MSG <- c(
    "Note that starting with BioC 3.7, the class attribute ",
    "of all GRangesList **instances** needs to be set to ",
    "\"CompressedGRangesList\". Please update this object ",
    "with 'updateObject(object, verbose=TRUE)' and ",
    "re-serialize it."
)

.valid.GenomicRangesList <- function(x)
{
    if (class(x) == "GRangesList")
        return(paste(OLD_GRANGESLIST_INSTANCE_MSG, collapse=""))
    .valid.GenomicRangesList.mcols(x)
}

setValidity2("GenomicRangesList", .valid.GenomicRangesList)


### TODO: Use this instead of GenomicRanges_OR_GRangesList in the
### RangedSummarizedExperiment class definition to specify the class
### of the 'rowRanges' slot.
setClassUnion("GenomicRanges_OR_GenomicRangesList",
    c("GenomicRanges", "GenomicRangesList")
)


setClass("SimpleGenomicRangesList",
    contains=c("GenomicRangesList", "SimpleRangesList"),
    representation("VIRTUAL")
)

setClass("CompressedGenomicRangesList",
    contains=c("GenomicRangesList", "CompressedRangesList"),
    representation("VIRTUAL")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###

### callNextMethod() searches for the "closest parent method" starting from
### the class used in the **signature** of the method from which it is called,
### NOT from 'class(x)'. This means that in the context of methods defined
### for GenomicRangesList objects, like the methods below, callNextMethod()
### would fail to find methods defined for CompressedList or SimpleList.
### So we will use .selectClosestMethod1(), which, unlike callNextMethod(),
### will find the "closest parent method" with respect to 'x'.
.selectClosestMethod1 <- function(f, x)
{
    if (is(x, "CompressedGenomicRangesList"))
        return(selectMethod(f, "CompressedRangesList"))
    if (is(x, "SimpleGenomicRangesList"))
        return(selectMethod(f, "SimpleRangesList"))
    stop(wmsg("'class(x)' (\"", class(x), "\") is not supported"))
}

setMethod("updateObject", "GenomicRangesList",
    function(object, ..., verbose=FALSE)
    {
        ## class attribute.
        if (class(object) == "GRangesList") {
            ## Starting with GenomicRanges 1.31.13, all GRangesList instances
            ## need to be replaced with CompressedGRangesList instances. Note
            ## that this is NOT a change of the internals (GRangesList
            ## instances have been using the CompressedList representation
            ## since the beginning), only a change of the class attribute.
            if (verbose)
                message("[updateObject] Settting class attribute of ",
                        "GRangesList instance\n",
                        "[updateObject] to \"CompressedGRangesList\" ... ",
                        appendLF=FALSE)
            class(object) <- class(new("CompressedGRangesList"))
            if (verbose)
                message("OK")
        } else {
            if (verbose)
                message("[updateObject] ", class(object), " object ",
                        "is current.\n",
                        "[updateObject] Nothing to update.")
        }

        if (is(object, "CompressedGenomicRangesList")) {
            ## unlistData slot.
            object@unlistData <- updateObject(object@unlistData,
                                              ..., verbose=verbose)
        }

        ## 'METHOD' will be the method for CompressedList objects or
        ## the method for Vector objects. (The former will update the
        ## 'partitioning' and 'elementMetadata' slots, while the latter
        ## will update the 'elementMetadata' slot only.)
        METHOD <- .selectClosestMethod1("updateObject", object)
        METHOD(object, ..., verbose=verbose)
    }
)

### Temporary hack to make 'length(x)' work on an old GRangesList instance.
setMethod("length", "GenomicRangesList",
    function(x)
    {
        if (class(x) == "GRangesList") {
            #warning(wmsg(OLD_GRANGESLIST_INSTANCE_MSG))
            x <- updateObject(x, check=FALSE)
        }
        METHOD <- .selectClosestMethod1("length", x)
        METHOD(x)
    }
)

### Temporary hack to make 'names(x)' work on an old GRangesList instance.
setMethod("names", "GenomicRangesList",
    function(x)
    {
        if (class(x) == "GRangesList") {
            #warning(wmsg(OLD_GRANGESLIST_INSTANCE_MSG))
            x <- updateObject(x, check=FALSE)
        }
        METHOD <- .selectClosestMethod1("names", x)
        METHOD(x)
    }
)

### Temporary hack to make 'names(x) <- value' work on an old GRangesList
### instance.
setReplaceMethod("names", "GenomicRangesList",
    function(x, value)
    {
        if (class(x) == "GRangesList") {
            #warning(wmsg(OLD_GRANGESLIST_INSTANCE_MSG))
            x <- updateObject(x, check=FALSE)
        }
        METHOD <- .selectClosestMethod1("names<-", x)
        METHOD(x, value)
    }
)

### Temporary hack to make 'extractROWS(x, i)' work on an old GRangesList
### instance.
setMethod("extractROWS", c("GenomicRangesList", "ANY"),
    function(x, i)
    {
        if (class(x) == "GRangesList") {
            #warning(wmsg(OLD_GRANGESLIST_INSTANCE_MSG))
            x <- updateObject(x, check=FALSE)
        }
        METHOD <- .selectClosestMethod1("extractROWS", x)
        METHOD(x, i)
    }
)

### Temporary hack to make 'getListElement(x, i)' work on an old GRangesList
### instance.
setMethod("getListElement", "GenomicRangesList",
    function(x, i, exact=TRUE)
    {
        if (class(x) == "GRangesList") {
            #warning(wmsg(OLD_GRANGESLIST_INSTANCE_MSG))
            x <- updateObject(x, check=FALSE)
        }
        METHOD <- .selectClosestMethod1("getListElement", x)
        METHOD(x, i, exact=exact)
    }
)

### Temporary hack to make 'unlist(x)' work on an old GRangesList instance.
setMethod("unlist", "GenomicRangesList",
    function(x, recursive=TRUE, use.names=TRUE)
    {
        if (class(x) == "GRangesList") {
            #warning(wmsg(OLD_GRANGESLIST_INSTANCE_MSG))
            x <- updateObject(x, check=FALSE)
        }
        if (is(x, "CompressedGenomicRangesList")) {
            if (!isTRUEorFALSE(use.names))
                stop("'use.names' must be TRUE or FALSE")
            ## Update any old GRanges instance (or other GenomicRanges
            ## derivative) stuck in the 'unlistData' slot.
            unlisted_x <- updateObject(x@unlistData, check=FALSE)
            if (use.names)
                unlisted_x <- S4Vectors:::set_unlisted_names(unlisted_x, x)
            return(unlisted_x)
        }
        METHOD <- .selectClosestMethod1("unlist", x)
        METHOD(x, recursive=recursive, use.names=use.names)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show()
###

### NOT exported but used in the GenomicAlignments package.
### FIXME: Seems to repeat a lot of code from IRanges:::show_IntegerRangesList!
show_GenomicRangesList <- function(x, with.header=TRUE)
{
    if (class(x) == "GRangesList") {
        warning(wmsg(OLD_GRANGESLIST_INSTANCE_MSG))
        x <- updateObject(x, check=FALSE)
    }
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

setMethod("seqinfo", "GenomicRangesList",
    combine_seqinfo_from_GenomicRanges_objects
)
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
            value <- make_zero_col_DFrame(length(x))
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
            value <- make_zero_col_DFrame(length(x@unlistData))
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
    if (!is(value, "Seqinfo"))
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    pruning.mode <- match.arg(pruning.mode)
    dangling_seqlevels <- GenomeInfoDb:::getDanglingSeqlevels(x,
                              new2old=new2old,
                              pruning.mode=pruning.mode,
                              seqlevels(value))
    if (length(dangling_seqlevels) != 0L) {
        ## Prune 'x'.
        idx <- !(seqnames(x) %in% dangling_seqlevels)
        if (pruning.mode == "coarse") {
            x <- x[all(idx)]  # "coarse" pruning
        } else {
            x <- x[idx]  # "fine" pruning
            if (pruning.mode == "tidy") {
                ## Remove list elements that became empty because of "fine"
                ## pruning.
                x <- x[any(idx) | elementNROWS(idx) == 0L]  # "tidy" pruning
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


