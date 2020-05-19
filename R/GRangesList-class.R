### =========================================================================
### GRangesList objects
### -------------------------------------------------------------------------
###


setClass("GRangesList",
    contains="GenomicRangesList",
    representation("VIRTUAL"),
    prototype(elementType="GRanges")
)

setClass("SimpleGRangesList",
    contains=c("GRangesList", "SimpleGenomicRangesList")
)

setClass("CompressedGRangesList",
    contains=c("GRangesList", "CompressedGenomicRangesList"),
    representation(unlistData="GRanges")
)

### TODO: Get rid of this! Used in RangedSummarizedExperiment class definition
### to specify the class of the 'rowRanges' slot and was originally introduced
### to support this use case. However, more packages use it these days e.g.
### DEFormats, GenomicFiles, ggbio, gmapR, HelloRanges, profileplyr, and maybe
### more... Everybody now should use GenomicRanges_OR_GenomicRangesList instead
### of GenomicRanges_OR_GRangesList.
setClassUnion("GenomicRanges_OR_GRangesList", c("GenomicRanges", "GRangesList"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

GRangesList <- function(..., compress=TRUE)
{
    if (!isTRUEorFALSE(compress))
        stop("'compress' must be TRUE or FALSE")
    objects <- list(...)
    if (length(objects) == 1L) {
        tmp <- objects[[1L]]
        if (is.list(tmp) || (is(tmp, "List") && !is(tmp, "GenomicRanges")))
            objects <- tmp
    }
    if (compress)
        suppressWarnings(as(objects, "CompressedGRangesList"))
    else
        as(objects, "SimpleGRangesList")
}

GenomicRangesList <- function(...)
{
    .Deprecated("GRangesList(..., compress=FALSE)")
    GRangesList(..., compress=FALSE)
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
    fragmentStarts <- normarg_list_of_integers(fragmentStarts, sep,
                                               "fragmentStarts")
    nfrag_per_feature <- elementNROWS(fragmentStarts)
    start <- unlist(fragmentStarts, recursive=FALSE, use.names=FALSE)

    fragmentEnds <- normarg_list_of_integers(fragmentEnds, sep,
                                             "fragmentEnds")
    nend_per_elt <- elementNROWS(fragmentEnds)
    if (length(nend_per_elt) != 0L) {
        if (length(nfrag_per_feature) == 0L)
            nfrag_per_feature <- nend_per_elt
        else if (!identical(nend_per_elt, nfrag_per_feature))
            stop("'fragmentStarts' and 'fragmentEnds' have ",
                 "incompatible \"shapes\"")
    }
    end <- unlist(fragmentEnds, recursive=FALSE, use.names=FALSE)

    fragmentWidths <- normarg_list_of_integers(fragmentWidths, sep,
                                               "fragmentWidths")
    nwidth_per_elt <- elementNROWS(fragmentWidths)
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
    partitioning <- PartitioningByEnd(cumsum(nfrag_per_feature), names=NULL)
    relist(unlistData, partitioning)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors.
###

setMethod("ranges", "CompressedGRangesList",
    function(x, use.names=TRUE, use.mcols=FALSE)
    {
        if (!isTRUEorFALSE(use.names))
            stop("'use.names' must be TRUE or FALSE")
        if (!isTRUEorFALSE(use.mcols))
            stop("'use.mcols' must be TRUE or FALSE")
        unlisted_x <- unlist(x, use.names=FALSE)
        unlisted_ans <- unlisted_x@ranges
        if (use.mcols)
            mcols(unlisted_ans) <- mcols(unlisted_x, use.names=FALSE)
        ans <- relist(unlisted_ans, x)
        if (!use.names)
            names(ans) <- NULL
        if (use.mcols)
            mcols(ans) <- mcols(x, use.names=FALSE)
        ans
    }
)

setReplaceMethod("ranges", "CompressedGRangesList",
    function(x, value)
    {
        if (!is(value, "IntegerRangesList") ||
            !identical(elementNROWS(x), elementNROWS(value)))
            stop("replacement 'value' is not an IntegerRangesList with the ",
                 "same elementNROWS as 'x'")
        ranges(x@unlistData) <- as(unlist(value, use.names = FALSE), "IRanges")
        x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion from list-like object to SimpleGRangesList
###

### Try to turn an arbitrary list-like object into an ordinary list of
### GRanges objects.
as_list_of_GRanges <- function(from)
{
    if (is(from, "GenomicRanges")) {
        if (!is(from, "GRanges"))
            from <- as(from, "GRanges", strict=FALSE)
        along_idx <- setNames(seq_along(from), names(from))
        names(from) <- NULL
        mcols(from) <- NULL
        lapply(along_idx, function(i) from[i])
    } else {
        lapply(from, as, "GRanges", strict=FALSE)
    }
}

### From ordinary list to SimpleGRangesList

.from_list_to_SimpleGRangesList <- function(from)
{
    from <- as_list_of_GRanges(from)
    S4Vectors:::new_SimpleList_from_list("SimpleGRangesList", from)
}

setAs("list", "SimpleGRangesList", .from_list_to_SimpleGRangesList)
setAs("list", "GRangesList", .from_list_to_SimpleGRangesList)

### From List derivative to SimpleGRangesList

.from_List_to_SimpleGRangesList <- function(from)
{
    S4Vectors:::new_SimpleList_from_list("SimpleGRangesList",
                               as_list_of_GRanges(from),
                               metadata=metadata(from),
                               mcols=mcols(from, use.names=FALSE))
}

setAs("List", "SimpleGRangesList", .from_List_to_SimpleGRangesList)

### Automatic coercion methods from SimpleList, GenomicRangesList, or
### SimpleGenomicRangesList to SimpleGRangesList silently return a broken
### object (unfortunately these dummy automatic coercion methods don't bother
### to validate the object they return). So we overwrite them.
setAs("SimpleList", "SimpleGRangesList",
      .from_List_to_SimpleGRangesList)
setAs("GenomicRangesList", "SimpleGRangesList",
      .from_List_to_SimpleGRangesList)
setAs("SimpleGenomicRangesList", "SimpleGRangesList",
      .from_List_to_SimpleGRangesList)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Intra-range methods
###

setMethod("trim", "GRangesList",
          function(x, use.names=TRUE)
          {
              ## Like seqinfo,GRangesList, assumes that there is a
              ## single Seqinfo for the entire object. Only guaranteed
              ## to be true in the compressed case.
              relist(trim(unlist(x, use.names=FALSE), use.names=use.names), x)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion from list-like object to CompressedGRangesList
###

### From ordinary list to CompressedGRangesList

.from_list_to_CompressedGRangesList <- function(from)
{
    from <- as_list_of_GRanges(from)
    IRanges:::new_CompressedList_from_list("CompressedGRangesList", from)
}

setAs("list", "CompressedGRangesList", .from_list_to_CompressedGRangesList)

### From List derivative to CompressedGRangesList

.from_List_to_CompressedGRangesList <- function(from)
{
    IRanges:::new_CompressedList_from_list("CompressedGRangesList",
                                 as_list_of_GRanges(from),
                                 metadata=metadata(from),
                                 mcols=mcols(from, use.names=FALSE))
}

### GenomicRanges objects are List objects so this case is already covered
### by the .from_List_to_CompressedGRangesList() helper above. However, we
### can implement it much more efficiently.
.from_GenomicRanges_to_CompressedGRangesList <- function(from)
{
    if (!is(from, "GRanges"))
        from <- as(from, "GRanges", strict=FALSE)
    ans_partitioning <- PartitioningByEnd(seq_along(from), names=names(from))
    names(from) <- NULL
    ans_mcols <- mcols(from, use.names=FALSE)
    mcols(from) <- NULL
    ans <- relist(from, ans_partitioning)
    mcols(ans) <- ans_mcols
    ans
}

setAs("List", "CompressedGRangesList",
      .from_List_to_CompressedGRangesList)

setAs("GenomicRanges", "CompressedGRangesList",
      .from_GenomicRanges_to_CompressedGRangesList)

setAs("List", "GRangesList",
    function(from)
    {
        if (is(from, "CompressedList") || is(from, "GenomicRanges"))
            as(from, "CompressedGRangesList")
        else
            as(from, "SimpleGRangesList")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other coercions
###

setAs("CompressedGRangesList", "CompressedIRangesList",
    function(from) ranges(from, use.mcols=TRUE)
)
setAs("CompressedGRangesList", "IRangesList",
    function(from) ranges(from, use.mcols=TRUE)
)
setAs("CompressedGRangesList", "IntegerRangesList",
    function(from) ranges(from, use.mcols=TRUE)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Sorting.
###
### S3/S4 combo for sort.GRangesList
.sort.GRangesList <- function(x, decreasing=FALSE, ...)
{
    gr <- deconstructGRLintoGR(x)
    gr2 <- sort(gr, decreasing=decreasing, ...)
    reconstructGRLfromGR(gr2, x)
}
sort.GRangesList <- function(x, decreasing=FALSE, ...)
    .sort.GRangesList(x, decreasing=decreasing, ...)
setMethod("sort", "CompressedGRangesList", .sort.GRangesList)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting.
###
### TODO: We should have more general methods defined at the GenomicRangesList
### level instead of the methods below.
###

.sBracketSubsetGRList <- function(x, i, j, ..., drop)
{
    if (!missing(i)) {
        x <- callNextMethod(x = x, i = i)
    }
    if (!missing(j)) {
        if (!is.character(j))
            stop("'j' must be a character vector")
        withinLevel <- (j %in% colnames(x@unlistData@elementMetadata))
        if (any(withinLevel) && !all(withinLevel))
            stop("'j' cannot mix between and within metadata column names")
        if (any(withinLevel)) {
            mcols(x, level="within") <-
              mcols(x, use.names=FALSE, level="within")[, j, drop=FALSE]
        } else {
            mcols(x) <- mcols(x, use.names=FALSE)[, j, drop=FALSE]
        }
    }
    x
}
setMethod("[", "CompressedGRangesList", .sBracketSubsetGRList)

.sBracketReplaceGRList <- function(x, i, j, ..., value)
{
    if (!is(value, class(x)[1]))
        stop(paste0("replacement value must be a ", class(x)[1], " object"))
    if (!missing(i))
        i <- extractROWS(setNames(seq_along(x), names(x)), i)
    if (!missing(j)) {
        if (!is.character(j))
            stop("'j' must be a character vector")
        withinLevel <- (j %in% colnames(x@unlistData@elementMetadata))
        if (any(withinLevel) && !all(withinLevel))
            stop("'j' cannot mix between and within metadata column names")
        if (missing(i)) {
            if (any(withinLevel)) {
                mcols(x, level="within")[, j] <-
                  mcols(x, use.names=FALSE, level="within")
            } else {
                mcols(x)[, j] <- mcols(x, use.names=FALSE)
            }
        } else {
            if (any(withinLevel)) {
                mcols(x, level="within")[i, j] <-
                        mcols(x, use.names=FALSE, level="within")
            } else {
                mcols(x)[i, j] <- mcols(x, use.names=FALSE)
            }
        }
    }
    callNextMethod(x = x, i = i, value = value)
}
setReplaceMethod("[", "CompressedGRangesList", .sBracketReplaceGRList)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Going from GRanges to GRangesList with extractList() and family.
###

setMethod("relistToClass", "GRanges", function(x) "CompressedGRangesList")

