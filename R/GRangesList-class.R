### =========================================================================
### GRangesList objects
### -------------------------------------------------------------------------
###

setClass("GRangesList",
    contains="CompressedList",
    prototype=prototype(
        elementType="GRanges",
        unlistData=new("GRanges"),
        elementMetadata=DataFrame()
    )
)


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
    if (length(listData) == 0L) {
        unlistData <- GRanges()
    } else {
        if (length(listData) == 1L && is.list(listData[[1L]]))
            listData <- listData[[1L]]
        if (!all(sapply(listData, is, "GRanges")))
            stop("all elements in '...' must be GRanges objects")
        unlistData <- suppressWarnings(do.call("c", unname(listData)))
    }
    end <- cumsum(elementLengths(unname(listData)))
    ans <- IRanges:::newCompressedList("GRangesList",
               unlistData,
               end = end, NAMES = names(listData),
               elementMetadata = new("DataFrame", nrows = length(listData)))
    validObject(ans)
    ans
}

setMethod("updateObject", "GRangesList",
    function(object, ..., verbose=FALSE)
    {
        if (verbose)
            message("updateObject(object = 'GRangesList')")
        if (!is(try(object@unlistData@seqinfo, silent=TRUE), "try-error"))
            return(object)
        object@unlistData <- updateObject(object@unlistData)
        object
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###

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

.GRangesListAsCompressedIRangesList <- function(from)
{
    ans_ranges <- from@unlistData@ranges
    ans_ranges@elementMetadata <- from@unlistData@elementMetadata
    new("CompressedIRangesList",
        unlistData=ans_ranges,
        partitioning=from@partitioning,
        elementMetadata=from@elementMetadata)
}

setAs("GRangesList", "CompressedIRangesList",
    .GRangesListAsCompressedIRangesList
)

setAs("GRangesList", "IRangesList",
    .GRangesListAsCompressedIRangesList
)

setAs("RangedDataList", "GRangesList",
      function(from) GRangesList(lapply(from, as, "GRanges")))


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

setMethod("seqinfo", "GRangesList", function(x) seqinfo(x@unlistData))

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
        seqnames(x@unlistData) <- value
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
        ranges(x@unlistData) <- as(unlist(value, use.names = FALSE), "IRanges")
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
        strand(x@unlistData) <- value
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
            x@elementMetadata <- value
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
            elementMetadata(x@unlistData) <- value
        }
        x
    }
)

setReplaceMethod("seqinfo", "GRangesList",
    function(x, value)
    {
        seqinfo(x@unlistData) <- value
        x
    }
)

setReplaceMethod("seqlevels", "GRangesList",
    function(x, value)
    {
        seqlevels(x@unlistData) <- value
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
    function(x, check = TRUE, value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        start(ranges(x@unlistData), check = check) <- value
        x
    }
)

setReplaceMethod("end", "GRangesList",
    function(x, check = TRUE, value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        end(ranges(x@unlistData), check = check) <- value
        x
    }
)

setReplaceMethod("width", "GRangesList",
    function(x, check = TRUE, value)
    {
        if (!is(value, "IntegerList") ||
            !identical(elementLengths(x), elementLengths(value)))
            stop("replacement 'value' is not an IntegerList with the same ",
                 "elementLengths as 'x'")
        value <- unlist(value, use.names = FALSE)
        width(ranges(x@unlistData), check = check) <- value
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
        ranges(x@unlistData) <-
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
### Sequence methods.
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

setReplaceMethod("[[", "GRangesList",
    function(x, i, j, ..., value)
    {
        nameValue <- if (is.character(i)) i else ""
        i <- IRanges:::normargSubset2_iOnly(x, i, j, ...,
                 .conditionPrefix="[[<-,GRangesList-method: ")
        len <- length(x)
        if (i > len) {
            value <- list(value)
            if (nzchar(nameValue))
                names(value) <- nameValue
            x <- c(x, do.call(getFunction(class(x)), value))
        } else {
            x <- callNextMethod(x, i, ..., value=value)
        }
        validObject(x)
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
            cat("<0 elements>\n\n")
        } else if ((k == 1L) || ((k <= 3L) && (N <= 20L))) {
            nms <- names(object)
            defnms <- paste("[[", seq_len(k), "]]", sep="")
            if (is.null(nms)) {
                nms <- defnms
            } else {
                empty <- nchar(nms) == 0L
                nms[empty] <- defnms[empty]
                nms[!empty] <- paste("$", nms[!empty], sep="")
            }
            for (i in seq_len(k)) {
                cat(nms[i], "\n")
                showGenomicRanges(object[[i]])
                cat("\n")
            }
        } else {
            sketch <- function(x) c(head(x, 3), "...", tail(x, 3))
            if (k >= 3 && cumsumN[3L] <= 20)
                showK <- 3
            else if (k >= 2 && cumsumN[2L] <= 20)
                showK <- 2
            else
                showK <- 1
            diffK <- k - showK
            nms <- names(object)[seq_len(showK)]
            defnms <- paste("[[", seq_len(showK), "]]", sep="")
            if (is.null(nms)) {
                nms <- defnms
            } else {
                empty <- nchar(nms) == 0L
                nms[empty] <- defnms[empty]
                nms[!empty] <- paste("$", nms[!empty], sep="")
            }
            for (i in seq_len(showK)) {
                cat(nms[i], "\n")
                showGenomicRanges(object[[i]])
                cat("\n")
            }
            if (diffK > 0) {
                cat("...\n<", k - showK,
                    ifelse(diffK == 1, " more element>\n\n", " more elements>\n\n"),
                    sep="")
            }
        }
        .showSeqlengths(object)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Deconstruction/reconstruction of a GRangesList into/from a GRanges
### object.
###
### For internal use only (not exported).
###

### Unlist GRangesList object 'x' into a GRanges object but the differences
### with the "unlist" method for GRangesList objects are:
###   - The sequence names of the returned GRanges object are modified by
###     embedding the "grouping by top-level element" information in them.
###   - The seqinfo is modified accordingly.
deconstructGRLintoGR <- function(x, expand.levels=FALSE)
{
    ans <- x@unlistData
    f1 <- rep.int(seq_len(length(x)), elementLengths(x))
    f2 <- as.integer(seqnames(ans))
    f12 <- paste(f1, f2, sep="|")

    ## Compute 'ans_seqinfo'.
    if (expand.levels) {
        x_nlev <- length(seqlevels(x))
        i1 <- rep(seq_len(length(x)), each=x_nlev)
        i2 <- rep.int(seq_len(x_nlev), length(x))
    } else {
        oo <- IRanges:::orderTwoIntegers(f1, f2)
        of1 <- f1[oo]
        of2 <- f2[oo]
        ## TODO: Add "presorted" method to IRanges:::duplicatedTwoIntegers()
        ## for when the 2 input vectors are already sorted.
        notdups <- !IRanges:::duplicatedTwoIntegers(of1, of2)
        i1 <- of1[notdups]
        i2 <- of2[notdups]
    }
    x_seqinfo <- seqinfo(x)
    ans_seqlevels <- paste(i1, i2, sep="|")
    ans_seqlengths <- unname(seqlengths(x_seqinfo))[i2]
    ans_isCircular <- unname(isCircular(x_seqinfo))[i2]
    ans_seqinfo <- Seqinfo(ans_seqlevels, ans_seqlengths, ans_isCircular)

    ## The 2 following modifications must be seen as a single atomic
    ## operation since doing the 1st without doing the 2nd would leave 'ans'
    ## in a broken state.
    ans@seqnames <- Rle(factor(f12, ans_seqlevels))
    ans@seqinfo <- ans_seqinfo
    ans
}

### The "inverse" transform of deconstructGRLintoGR().
### More precisely, reconstructGRLfromGR() transforms GRanges object 'gr'
### with sequence names in the "f1|f2" format (as produced by
### deconstructGRLintoGR() above) back into a GRangesList object with the
### same length & names & elementMetadata & seqinfo as 'x'.
### The fundamental property of this deconstruction/reconstruction mechanism
### is that, for any GRangesList object 'x':
###
###   reconstructGRLfromGR(deconstructGRLintoGR(x), x) is identical to x
###
reconstructGRLfromGR <- function(gr, x)
{
    snames <- strsplit(as.character(seqnames(gr)), "|", fixed=TRUE)
    m12 <- matrix(as.integer(unlist(snames)), ncol=2, byrow=TRUE)

    ## Restore the real sequence names.
    f2 <- m12[ , 2L]
    x_seqlevels <- seqlevels(x)
    ## The 2 following modifications must be seen as a single atomic
    ## operation since doing the 1st without doing the 2nd would leave 'ans'
    ## in a broken state.
    gr@seqnames <- Rle(factor(x_seqlevels[f2], x_seqlevels))
    gr@seqinfo <- seqinfo(x)

    ## Split.
    f1 <- m12[ , 1L]
    ans <- split(gr, factor(f1, levels=seq_len(length(x))))
    names(ans) <- names(x)
    elementMetadata(ans) <- elementMetadata(x)
    ans
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The "range" and "reduce" methods.
###
### For a GRangesList object 'x', 'range(x)' and 'reduce(x)' are equivalent
### to 'endoapply(x, range)' and 'endoapply(x, reduce)', respectively.
### This makes them isomorphisms, that is, they are endomorphisms (i.e. they
### preserve the class of 'x') who also preserve the length & names &
### elementMetadata of 'x'. In addition, the seqinfo is preserved too.
###
### However, using endoapply() for the implementation would be too
### inefficient. The fast implementation below takes advantage of the
### fact that we already have fast "range" and "reduce" methods for GRanges
### objects. Depending on the size of 'x', the implementation below is 50x or
### 1000x faster (or more) than the implementation using endoapply().

setMethod("range", "GRangesList",
    function(x, ..., na.rm=FALSE)
    {
        if (length(list(...)) != 0L)
            stop("\"range\" method for GRangesList objects only ",
                 "takes a single object")
        gr <- deconstructGRLintoGR(x)
        ## "range" method for GRanges objects is fast.
        gr <- range(gr)
        reconstructGRLfromGR(gr, x)
    }
)

setMethod("reduce", "GRangesList",
    function(x, drop.empty.ranges=FALSE, min.gapwidth=1L,
             with.inframe.attrib=FALSE)
    {
        if (!identical(with.inframe.attrib, FALSE)) 
            stop("'with.inframe.attrib' argument is not supported ", 
                 "when reducing a GRangesList object")
        gr <- deconstructGRLintoGR(x)
        ## "reduce" method for GRanges objects is fast.
        gr <- reduce(gr, drop.empty.ranges=drop.empty.ranges,
                         min.gapwidth=min.gapwidth)
        reconstructGRLfromGR(gr, x)
    }
)

