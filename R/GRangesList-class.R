### =========================================================================
### GRangesList objects
### -------------------------------------------------------------------------
###
### Class definition

setClass("GRangesList", contains = "CompressedList",
         prototype = prototype(elementType = "GRanges",
                               unlistData = new("GRanges"),
                               elementMetadata = DataFrame()))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity.
###

.valid.GRangesList.elementMetadata <- function(x)
{
    msg <- NULL
    if (nrow(x@elementMetadata) != length(x))
        msg <- "slot 'elementMetadata' has an incorrect number of rows"
    if (any(c("seqnames", "ranges", "strand", "start", "end", "width",
              "element") %in% colnames(x@elementMetadata)))
        msg <-
          c(msg,
            paste("slot 'elementMetadata' cannot use \"seqnames\", \"ranges\",",
                  "\"strand\", \"start\", \"end\", \"width\", or \"element\"",
                  "as column names"))
    if (any(colnames(x@elementMetadata) %in%
            colnames(x@unlistData@elementMetadata)))
        msg <-
          c(msg,
            paste("slot 'elementMetadata' cannot use the same names for",
                  "columns as unlisted GRanges elementMetadata"))
    if (!is.null(rownames(x@elementMetadata)))
        msg <- c(msg, "slot 'elementMetadata' cannot contain row names")
    msg
}

.valid.GRangesList <- function(x)
{
    c(.valid.GRangesList.elementMetadata(x))
}

setValidity2("GRangesList", .valid.GRangesList)


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
    ans <-
      IRanges:::newCompressedList("GRangesList", listData,
                                  elementMetadata =
                                  new("DataFrame", nrows = length(listData)))
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
            element <- rep.int(seq_len(length(x)), elementLengths(x))
        else
            element <- rep.int(names(x), elementLengths(x))
        data.frame(element = element,
                   as.data.frame(unlist(x, use.names = FALSE),
                                 row.names = row.names),
                   stringsAsFactors = FALSE)
    }
)

setMethod("unlist", "GRangesList",
    function(x, recursive = TRUE, use.names = TRUE)
    {
        if (use.names && is.null(names(x@unlistData)))
            names(x@unlistData) <- seq_len(length(x@unlistData))
        callNextMethod()
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

setMethod("elementMetadata", "GRangesList",
    function(x, level = c("between", "within"), ...)
    {
        level <- match.arg(level)
        if (level == "between") {
            ans <- x@elementMetadata
            if (!is.null(names(x)))
                rownames(ans) <- names(x)
            ans
        } else {
            elementMetadata <- x@unlistData@elementMetadata
            nms <- names(x@unlistData)
            if (!is.null(nms))
                rownames(elementMetadata) <- nms
            ans <-
              new2("CompressedSplitDataFrameList", unlistData = elementMetadata,
                   partitioning = x@partitioning, check=FALSE)
        }
        ans
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
            value <- Rle(factor(value))
        else if (!is.factor(runValue(value)))
            runValue(value) <- factor(runValue(value))
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

setReplaceMethod("elementMetadata", "GRangesList",
    function(x, level = c("between", "within"), ..., value) 
    {
        level <- match.arg(level)
        if (level == "between") {
            if (is.null(value))
                value <- new("DataFrame", nrows = length(x))
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
            x@elementData <- value
        } else {
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
            x@unlistData@elementMetadata <- value
        }
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
        if (!missing(i)) {
            x <- callNextMethod(x = x, i = i)
        }
        if (!missing(j)) {
            if (!is.character(j))
                stop("'j' must be a character vector")
            withinLevel <- (j %in% colnames(x@unlistData@elementMetadata))
            if (any(withinLevel) && !all(withinLevel))
                stop("'j' cannot mix between and within column names")
            if (any(withinLevel)) {
                elementMetadata(x, level="within") <-
                  elementMetadata(x, level="within")[, j, drop=FALSE]
            } else {
                elementMetadata(x) <- elementMetadata(x)[, j, drop=FALSE]
            }
        }
        x
    }
)

setReplaceMethod("[", "GRangesList",
    function(x, i, j, ..., value)
    {
        if (!is(value, "GRangesList"))
            stop("replacement value must be a GRangesList object")
        if (!missing(i)) {
            iInfo <- IRanges:::.bracket.Index(i, length(x), names(x))
            if (!is.null(iInfo[["msg"]]))
                stop(iInfo[["msg"]])
            i <- iInfo[["idx"]]
        }
        if (!missing(j)) {
            if (!is.character(j))
                stop("'j' must be a character vector")
            withinLevel <- (j %in% colnames(x@unlistData@elementMetadata))
            if (any(withinLevel) && !all(withinLevel))
                stop("'j' cannot mix between and within column names")
            if (missing(i)) {
                if (any(withinLevel)) {
                    elementMetadata(x, level="within")[, j] <-
                      elementMetadata(x, level="within")
                } else {
                    elementMetadata(x)[, j] <- elementMetadata(x)
                }
            } else {
                if (any(withinLevel)) {
                    elementMetadata(x, level="within")[i, j] <-
                            elementMetadata(x, level="within")
                } else {
                    elementMetadata(x)[i, j] <- elementMetadata(x)
                }
            }
        }
        callNextMethod(x = x, i = i, value = value)
    }
)

setMethod("ncol", "GRangesList", function(x) ncol(x@elementMetadata))

setMethod("colnames", "GRangesList",
    function(x, do.NULL = TRUE, prefix = "col") 
        colnames(x@elementMetadata, do.NULL = do.NULL, prefix = prefix))
setReplaceMethod("colnames", "GRangesList",
    function(x, value)
    {
        colnames(x@elementMetadata) <- value
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
