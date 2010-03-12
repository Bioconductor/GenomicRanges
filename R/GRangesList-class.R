### =========================================================================
### GRangesList objects
### -------------------------------------------------------------------------
###
### Class definition

setClass("GRangesList", contains = "CompressedList",
         prototype = prototype(elementType = "GRanges",
                               unlistData = new("GRanges")))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###

GRangesList <- function(...)
{
    listData <- list(...)
    if (length(listData) == 1 && is.list(listData[[1L]]))
        listData <- listData[[1L]]
    if (!all(sapply(listData, is, "GRanges")))
        stop("all elements in '...' must be GRanges objects")
    ans <- IRanges:::newCompressedList("GRangesList", listData)
    validObject(ans)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

setAs("RangedDataList", "GRangesList",
      function(from) GRangesList(lapply(from, as, "GRanges")))

setMethod("as.data.frame", "GRangesList",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        if (missing(row.names))
            row.names <- names(x@unlistData)
        if (is.null(names(x)))
            feature <- rep(seq_len(length(x)), elementLengths(x))
        else
            feature <- rep(names(x), elementLengths(x))
        data.frame(feature = feature,
                   as.data.frame(unlist(x, use.names = FALSE),
                                 row.names = row.names),
                   stringsAsFactors = FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods.
###

setMethod("seqnames", "GRangesList",
    function(x)
        new2("CompressedRleList",
             unlistData = x@unlistData@seqnames, partitioning = x@partitioning,
             check=FALSE))

setMethod("ranges", "GRangesList",
    function(x, ...)
        new2("CompressedIRangesList",
             unlistData = x@unlistData@ranges, partitioning = x@partitioning,
             check=FALSE))

setMethod("strand", "GRangesList",
    function(x)
        new2("CompressedRleList",
             unlistData = x@unlistData@strand, partitioning = x@partitioning,
             check=FALSE))

setMethod("values", "GRangesList",
    function(x, ...)
    {
        values <- x@unlistData@values
        nms <- names(x@unlistData)
        if (!is.null(nms))
            rownames(values) <- nms
        new2("CompressedSplitDataFrameList",
             unlistData = values, partitioning = x@partitioning, check=FALSE)
    }
)

setReplaceMethod("seqnames", "GRangesList",
    function(x, value) 
    {
        if (!is(value, "AtomicList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an AtomicList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        if (!is(value, "Rle"))
            value <- Rle(as.character(value))
        else if (!is.character(runValue(value)))
            runValue(value) <- as.character(runValue)
        x@unlistData@seqnames <- value
        x
    }
)

setReplaceMethod("ranges", "GRangesList",
    function(x, value) 
    {
        if (!is(value, "RangesList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not a RangesList with the same ",
                 "elementLengths as 'x'")
        x@unlistData@ranges <- as(unlist(value, use.names = FALSE), "IRanges")
        x
    }
)

setReplaceMethod("strand", "GRangesList",
    function(x, value) 
    {
        if (!is(value, "AtomicList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an AtomicList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        if (!is(value, "Rle"))
            value <- Rle(strand(value))
        else if (!is.factor(runValue(value)) ||
                 !identical(levels(runValue(value)), levels(strand())))
            runValue(value) <- strand(runValue(value))
        x@unlistData@strand <- value
        x
    }
)

setReplaceMethod("values", "GRangesList",
    function(x, value) 
    {
        if (is.null(value)) {
            value <- new("DataFrame", nrows = length(x@unlistData))
        } else {
            if (!is(value, "SplitDataFrameList") ||
                !identical(elementLengths(x), elementLengths(value))) {
                stop("replacement 'value' is not a SplitDataFrameList with ",
                     "the same elementLengths as 'x'")
             }
             value <- unlist(value, use.names = FALSE)
        }
        x@unlistData@values <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### RangesList methods.
###

setMethod("start", "GRangesList",
    function(x, ...)
        new2("CompressedIntegerList",
             unlistData = start(x@unlistData@ranges),
             partitioning = x@partitioning, check=FALSE))

setMethod("end", "GRangesList",
    function(x, ...)
        new2("CompressedIntegerList",
             unlistData = end(x@unlistData@ranges),
             partitioning = x@partitioning, check=FALSE))

setMethod("width", "GRangesList",
    function(x)
        new2("CompressedIntegerList",
             unlistData = width(x@unlistData@ranges),
             partitioning = x@partitioning, check=FALSE))

setReplaceMethod("start", "GRangesList",
    function(x, value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        start(x@unlistData@ranges) <- value
        x
    }
)

setReplaceMethod("end", "GRangesList",
    function(x, value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        end(x@unlistData@ranges) <- value
        x
    }
)

setReplaceMethod("width", "GRangesList",
    function(x, value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        width(x@unlistData@ranges) <- value
        x
    }
)

setMethod("shift", "GRangesList",
    function(x, shift, use.names=TRUE)
    {
        if (is(shift, "IntegerList")) {
            if (length(shift) != length(x) ||
                any(elementLengths(shift) != elementLengths(x))) {
                stop("IntegerList 'shift' not of same dimension as 'x'")
            }
            shift <- unlist(shift, use.names=FALSE)
        }
        x@unlistData@ranges <-
          shift(x@unlistData@ranges, shift, use.names=use.names)
        x
    }
)

setMethod("coverage", "GRangesList",
    function(x, shift = list(0L), width = list(NULL), weight = list(1L))
    {
        callGeneric(x@unlistData, shift = shift, width = width, weight = weight)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SplitDataFrameList methods.
###

setMethod("[", "GRangesList",
    function(x, i, j, ..., drop)
    {
        if (!missing(i))
            x <- callNextMethod(x = x, i = i)
        if (!missing(j))
            values(x) <- values(x)[, j, drop=FALSE]
        x
    }
)

setReplaceMethod("[", "GRangesList",
    function(x, i, j, ..., value)
    {
        if (!is(value, "GRangesList"))
            stop("replacement value must be a GRangesList object")
        if (!missing(j)) {
            subvalues <- values(x)[i,,drop=FALSE]
            subvalues[,j] <- values(value)
            values(value) <- subvalues
        }
        callNextMethod(x = x, i = i, value = value)
    }
)

setMethod("ncol", "GRangesList", function(x) ncol(x@unlistData@values))

setMethod("colnames", "GRangesList",
    function(x, do.NULL = TRUE, prefix = "col") 
        colnames(x@unlistData@values, do.NULL = do.NULL, prefix = prefix))
setReplaceMethod("colnames", "GRangesList",
    function(x, value)
    {
        colnames(x@unlistData@values) <- value
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show method.
###

setMethod("show", "GRangesList",
    function(object)
    {
        k <- length(object)
        cumsumN <- cumsum(elementLengths(object))
        N <- tail(cumsumN, 1)
        cat(class(object), " of length ", k, "\n", sep = "")
        if (k == 0L) {
            cat("<0 elements>\n")
        } else if ((k == 1L) || (N <= 20L)) {
            show(as.list(object))
        } else {
            sketch <- function(x) c(head(x, 3), "...", tail(x, 3))
            if (k >= 3 && cumsumN[3L] <= 20)
                showK <- 3
            else if (k >= 2 && cumsumN[2L] <= 20)
                showK <- 2
            else
                showK <- 1
            diffK <- k - showK
            show(as.list(object[seq_len(showK)]))
            if (diffK > 0)
                cat("...\n<", k - showK,
                    ifelse(diffK == 1, " more element>\n", " more elements>\n"),
                    sep="")
        }
    }
)
