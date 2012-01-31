### =========================================================================
### GRangesList objects
### -------------------------------------------------------------------------
###

setClass("GRangesList",
    contains=c("CompressedList", "GenomicRangesList"),
    representation(
        unlistData="GRanges",
        elementMetadata="DataFrame"
    ),
    prototype(
        elementType="GRanges"
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
### Constructors.
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

### TODO: Use this in GenomicFeatures::transcriptLocs2refLocs() and remove
### GenomicFeatures:::.normargExonStartsOrEnds().
.normargListOfIntegers <- function(arg, sep, argname) 
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

### Typically, the field values will come from a file that needs to be loaded
### into a data.frame first.
makeGRangesListFromFeatureFragments <- function(seqnames=Rle(factor()),
                                                fragmentStarts=list(),
                                                fragmentEnds=list(),
                                                fragmentWidths=list(),
                                                strand=character(0),
                                                sep=",")
{
    fragmentStarts <- .normargListOfIntegers(fragmentStarts, sep,
                                             "fragmentStarts")
    nfrag_per_feature <- elementLengths(fragmentStarts)
    start <- unlist(fragmentStarts, recursive=FALSE, use.names=FALSE)

    fragmentEnds <- .normargListOfIntegers(fragmentEnds, sep,
                                           "fragmentEnds")
    nend_per_elt <- elementLengths(fragmentEnds)
    if (length(nend_per_elt) != 0L) {
        if (length(nfrag_per_feature) == 0L)
            nfrag_per_feature <- nend_per_elt
        else if (!identical(nend_per_elt, nfrag_per_feature))
            stop("'fragmentStarts' and 'fragmentEnds' have ",
                 "incompatible \"shapes\"")
    }
    end <- unlist(fragmentEnds, recursive=FALSE, use.names=FALSE)

    fragmentWidths <- .normargListOfIntegers(fragmentWidths, sep,
                                             "fragmentWidths")
    nwidth_per_elt <- elementLengths(fragmentWidths)
    if (length(nwidth_per_elt) != 0L) {
        if (length(nfrag_per_feature) == 0L)
            nfrag_per_feature <- nwidth_per_elt
        else if (!identical(nwidth_per_elt, nfrag_per_feature))
            stop("\"shape\" of 'fragmentWidths' is incompatible ",
                 "with \"shape\" of 'fragmentStarts' or 'fragmentEnds'")
    }
    width <- unlist(fragmentWidths, recursive=FALSE, use.names=FALSE)

    ranges <- IRanges(start=start, end=end, width=width)
    nfrag <- sum(nfrag_per_feature)
    if (nfrag != length(ranges))
        stop("GenomicRanges internal error in makeGRangesListFromFields(): ",
             "nfrag != length(ranges). This should never happen. ",
             "Please report.")
    if (nfrag == 0L) {
        ## Cannot blindly subset by FALSE because it doesn't work on a
        ## zero-length Rle.
        if (length(seqnames) != 0L)
            seqnames <- seqnames[FALSE]
        if (length(strand) != 0L)
            strand <- strand[FALSE]
    } else {
        if (length(seqnames) != length(nfrag_per_feature) ||
            length(strand) != length(nfrag_per_feature))
            stop("length of 'seqnames' and/or 'strand' is incompatible ",
                 "with fragmentStarts/Ends/Widths")
        seqnames <- rep.int(seqnames, nfrag_per_feature)
        strand <- rep.int(strand, nfrag_per_feature)
    }
    unlistData <- GRanges(seqnames=seqnames, ranges=ranges, strand=strand)
    ans <- IRanges:::newCompressedList("GRangesList",
             unlistData, end=cumsum(nfrag_per_feature), NAMES=NULL,
             elementMetadata=new("DataFrame", nrows=length(nfrag_per_feature)))
    validObject(ans)
    ans
}

setMethod("updateObject", "GRangesList",
    function(object, ..., verbose=FALSE)
    {
        if (verbose)
            message("updateObject(object = 'GRangesList')")
        if (is(try(validObject(object@unlistData, complete=TRUE), silent=TRUE),
               "try-error")) {
            object@unlistData <- updateObject(object@unlistData)
            return(object)
        }
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
    function(x, row.names=FALSE, level=c("between", "within"), ...)
    {
        if (!isTRUEorFALSE(row.names))
            stop("'row.names' must be TRUE or FALSE")
        level <- match.arg(level)
        if (level == "between") {
            ans <- x@elementMetadata
            if (row.names)
                rownames(ans) <- names(x)
        } else {
            elementMetadata <- x@unlistData@elementMetadata
            if (row.names)
                rownames(elementMetadata) <- names(x@unlistData)
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
    function(x, new2old=NULL, force=FALSE, value)
    {
        if (!is(value, "Seqinfo"))
            stop("the supplied 'seqinfo' must be a Seqinfo object")
        unlisted <- x@unlistData
        dangling_seqlevels <- getDanglingSeqlevels(unlisted,
                                  new2old=new2old, force=force,
                                  seqlevels(value))
        if (length(dangling_seqlevels) != 0L) {
            dropme <- which(seqnames(unlisted) %in% dangling_seqlevels)
            x <- x[-unique(togroup(x, j=dropme))]
        }
        seqinfo(x@unlistData, new2old=new2old) <- value
        x
    }
)

setMethod("score", "GRangesList", function(x) {
  values(x)$score
})

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
    function(x, shift=0L, use.names=TRUE)
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

setMethod("isDisjoint", "GRangesList",
    function(x, ignore.strand = FALSE)
    {
        gr <- GenomicRanges:::deconstructGRLintoGR(x)

        if (ignore.strand) 
            xIRangesList <- split(unname(ranges(gr)), paste(seqnames(gr),
                           Rle(factor(rep("+", length(gr)))), sep = "\r"))
        else 
            xIRangesList <- split(unname(ranges(gr)),
                                  paste(seqnames(gr), strand(gr), sep = "\r"))
         
        ansIRanges <- isDisjoint(xIRangesList)
        splitListNames <- strsplit(names(ansIRanges), split="\r")
        snames <- strsplit(unlist(lapply(splitListNames, "[[", 1L)), "|",
                           fixed=TRUE)
        m12 <- matrix(as.integer(unlist(snames)), ncol=2, byrow=TRUE)
        
        ansIRangesList <- split(ansIRanges,
                                factor(m12[, 1L], levels=seq_len(length(x))))
        ans <-  unlist(lapply(ansIRangesList, all))
        names(ans) <- names(x)
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Vector methods.
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
        cat(class(object), " of length ", k, ":\n", sep = "")
        with.classinfo <- TRUE
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
                showGenomicRanges(object[[i]], margin="  ",
                                  with.classinfo=with.classinfo)
                if (with.classinfo)
                    with.classinfo <- FALSE
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
                showGenomicRanges(object[[i]], margin="  ",
                                  with.classinfo=with.classinfo)
                if (with.classinfo)
                    with.classinfo <- FALSE
                cat("\n")
            }
            if (diffK > 0) {
                cat("...\n<", k - showK,
                    ifelse(diffK == 1, " more element>\n", " more elements>\n"),
                    sep="")
            }
        }
        cat("---\n")
        showSeqlengths(object)
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
        oo <- IRanges:::orderIntegerPairs(f1, f2)
        of1 <- f1[oo]
        of2 <- f2[oo]
        ## TODO: Add "presorted" method to IRanges:::duplicatedIntegerPairs()
        ## for when the 2 input vectors are already sorted.
        notdups <- !IRanges:::duplicatedIntegerPairs(of1, of2)
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
    metadata(ans) <- metadata(x)
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

setMethod("restrict", "GRangesList",
    function(x, start = NA, end = NA, keep.all.ranges = FALSE, use.names = TRUE)
    {
        endoapply(x, restrict, start=start, end=end,keep.all.ranges=keep.all.ranges
               , use.names=use.names )

    })

setMethod("flank", "GRangesList",
    function(x, width, start=TRUE, both=FALSE, use.names=TRUE, ignore.strand=FALSE)
    {
        x@unlistData <- flank(x@unlistData, width=width, start=start, both=FALSE, use.names=TRUE,
                    ignore.strand=ignore.strand)
        x
    })


.listCumsum <- function(x) {
  x_unlisted <- unlist(x, use.names=FALSE)
  x_cumsum <- cumsum(x_unlisted)
  x_part <- PartitioningByWidth(elementLengths(x))
  x_cumsum - rep(x_cumsum[start(x_part)] - x_unlisted[start(x_part)],
                 width(x_part))
}

.listCumsumShifted <- function(x) {
  cs <- .listCumsum(x)
  shifted <- c(0L, head(cs, -1))
  shifted[start(PartitioningByWidth(elementLengths(x)))] <- 0L
  shifted
}

setMethod("map", c("GenomicRanges", "GRangesList"), function(from, to) {
  gr <- unlist(to, use.names=FALSE)
  ol <- findOverlaps(from, gr, type = "within")
  shits <- subjectHits(ol)
  qhits <- queryHits(ol)
  local <- ranges(from)[qhits]
  bounds <- ranges(gr)[shits]
  
  ## location wrt start of coding region 
  neg <- as.vector(strand(gr)[shits] == "-")
  local[!neg] <- shift(local[!neg], - start(bounds)[!neg])
  local[neg] <- IRanges(end(bounds)[neg] - end(local)[neg],
                        width = width(local)[neg])
  
  ## location wrt transcript 
  cumsums <- .listCumsumShifted(width(to))
  local <- shift(local, 1L + cumsums[shits])
  
  toInd <- rep(seq(length(to)), elementLengths(to))[shits]
  matching <- new("RangesMatching",
                  matchMatrix = cbind(query = qhits, subject = toInd),
                  queryLength = length(from),
                  subjectLength = length(to))
  new("RangesMapping", matching = matching, ranges = local,
      space = seqnames(from)[qhits])
})
