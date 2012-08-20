setClassUnion("GenomicRangesORGRangesList", c("GenomicRanges", "GRangesList"))

setClass("SummarizedExperiment",
    representation(
        exptData="SimpleList",                # overall description
        rowData="GenomicRangesORGRangesList", # rows and their annotations
        colData="DataFrame",                  # columns and their annotations
        assays="SimpleList"),                 # Data -- e.g., list of matricies
    prototype(rowData=GRanges()))

.valid.SummarizedExperiment <- function(x)
{
    msg <- NULL
    msg1 <-
        "\n   assays(<SummarizedExperiment>)[[%d]] must be a matrix or array" 
    msg2 <-
        "\n   dim(<SummarizedExperiment>) does not equal dim(assays(<SummarizedExperiment>)[[%d]])"

    xdim <- dim(x)
    assays <- assays(x)
    for (i in seq_len(length(assays))) {
        if (!is(assays[[i]], "matrix") & !is(assays[[i]], "array")) {
            msg <- c(msg, sprintf(msg1, i))
            next
        }
        edim <- dim(assays[[i]])[1:2]
        if (is.null(edim) || !all(xdim == edim))
            msg <- c(msg, sprintf(msg2, i))
    }

    msg
}

setValidity2("SummarizedExperiment", .valid.SummarizedExperiment)

setGeneric("SummarizedExperiment",
    function(assays, ...) standardGeneric("SummarizedExperiment"))

.GRangesList_assays <-
    function(assays)
{
    m <- assays[[1]]
    n <- nrow(m)
    names <- rownames(m)
    relist(GRanges(), PartitioningByEnd(integer(n), names=names))
}

setMethod(SummarizedExperiment, "SimpleList",
    function(assays, rowData=GRangesList(), colData=DataFrame(),
             exptData=SimpleList(), ...,
             verbose=FALSE)
{
    if (missing(rowData) && 0L != length(assays))
        rowData <- .GRangesList_assays(assays)
    if (missing(colData) && 0L != length(assays)) {
        nms <- colnames(assays[[1]])
        if (is.null(nms) && 0L != ncol(assays[[1]]))
            stop("SummarizedExperiment assay colnames must not be NULL")
        colData <- DataFrame(row.names=nms)
    }
    ## FIXME: warn if dimnames(assays) != list( verbose=TRUE
    if (!all(sapply(assays, function(x) is.null(dimnames(x)))))
        assays <- endoapply(assays, "dimnames<-", NULL)
    new("SummarizedExperiment", exptData=exptData,
        rowData=rowData, colData=colData, assays=assays, ...)
})

setMethod(SummarizedExperiment, "missing",
    function(assays, ...)
{
    SummarizedExperiment(SimpleList(), ...)
})

setMethod(SummarizedExperiment, "list",
    function(assays, ...)
{
    SummarizedExperiment(do.call(SimpleList, assays), ...)
})

setMethod(SummarizedExperiment, "matrix",
    function(assays, ...)
{
    SummarizedExperiment(SimpleList(assays), ...)
})

## FIXME: MTM thinks that generics do not belong with class definitions

setGeneric("exptData", function(x, ...) standardGeneric("exptData"))

setGeneric("exptData<-",
    function(x, ..., value) standardGeneric("exptData<-"))

## rowData, colData seem too vague, but from eSet derived classes wanted to
## call the rows / cols something different from 'features' or 'samples', so
## might as well avoid the issue
setGeneric("rowData", function(x, ...) standardGeneric("rowData"))

setGeneric("rowData<-",
    function(x, ..., value) standardGeneric("rowData<-"))

setGeneric("colData", function(x, ...) standardGeneric("colData"))

setGeneric("colData<-",
    function(x, ..., value) standardGeneric("colData<-"))

setGeneric("assays",
    function(x, ..., withDimnames=TRUE) standardGeneric("assays"),
    signature="x")

setGeneric("assays<-",
    function(x, ..., value) standardGeneric("assays<-"))

setGeneric("assay", function(x, i, ...) standardGeneric("assay"))

setGeneric("assay<-",
    function(x, i, ..., value) standardGeneric("assay<-"))

## Simple 'getters' / 'setters'

setMethod(exptData, "SummarizedExperiment",
    function(x, ...) slot(x, "exptData"))

setReplaceMethod("exptData", c("SummarizedExperiment", "SimpleList"),
    function(x, ..., value)
{
    initialize(x, ..., exptData=value)
})

setReplaceMethod("exptData", c("SummarizedExperiment", "list"),
    function(x, ..., value)
{
    initialize(x, ..., exptData=SimpleList(value))
})

setMethod(rowData, "SummarizedExperiment",
    function(x, ...) slot(x, "rowData"))

setReplaceMethod("rowData", c("SummarizedExperiment", "GenomicRanges"),
    function(x, ..., value)
{
    initialize(x, ..., rowData=value)
})

setReplaceMethod("rowData", c("SummarizedExperiment", "GRangesList"),
    function(x, ..., value)
{
    initialize(x, ..., rowData=value)
})

## also seqlevels, genome, seqlevels<-, genome<-
setMethod(seqinfo, "SummarizedExperiment",
    function(x)
{
    seqinfo(rowData(x))
})

setReplaceMethod("seqinfo", "SummarizedExperiment",
    function(x, new2old = NULL, force = FALSE, value)
{
    seqinfo(rowData(x), new2old=new2old, force=force) <- value
    x
})

setMethod(mcols, "SummarizedExperiment",
    function(x, use.names=FALSE, ...)
{
    mcols(rowData(x), use.names=use.names, ...)
})

setReplaceMethod("mcols", "SummarizedExperiment",
    function(x, ..., value)
{
    mcols(rowData(x), ...) <- value
    x
})

### mcols() is the recommended way for accessing the metadata columns.
### Use of values() or elementMetadata() is discouraged.

setMethod(elementMetadata, "SummarizedExperiment",
    function(x, use.names=FALSE, ...)
{
    elementMetadata(rowData(x), use.names=use.names, ...)
})

setReplaceMethod("elementMetadata", "SummarizedExperiment",
    function(x, ..., value)
{
    elementMetadata(rowData(x), ...) <- value
    x
})

setMethod(values, "SummarizedExperiment",
    function(x, ...)
{
    values(rowData(x), ...)
})

setReplaceMethod("values", "SummarizedExperiment",
    function(x, ..., value)
{
    values(rowData(x), ...) <- value
    x
})

setMethod(colData, "SummarizedExperiment",
    function(x, ...) slot(x, "colData"))

setReplaceMethod("colData", c("SummarizedExperiment", "DataFrame"),
    function(x, ..., value)
{
    initialize(x, ..., colData=value)
})

setMethod(assays, "SummarizedExperiment",
function(x, ..., withDimnames=TRUE) 
{
    if (withDimnames)
        endoapply(slot(x, "assays"), "dimnames<-", dimnames(x))
    else
        slot(x, "assays")
})

setReplaceMethod("assays", c("SummarizedExperiment", "SimpleList"),
    function(x, ..., value)
{
    initialize(x, ..., assays=value)
})

setReplaceMethod("assays", c("SummarizedExperiment", "list"),
    function(x, ..., value)
{
    value <- lapply(x, "dimnames<-", NULL)
    initialize(x, ..., assays=SimpleList(value))
})

## convenience for common use case 
setMethod(assay, c("SummarizedExperiment", "missing"),
    function(x, i, ...)
{
    assays <- assays(x, ...)
    if (0L == length(assays))
    {
        msg <- 'assay(<SummarizedExperiment>, i="missing", ...) length(assays(<SummarizedExperiment>)) is 0'
        stop(msg)
    }
    assays[[1]]
})

setMethod(assay, c("SummarizedExperiment", "numeric"),
    function(x, i, ...) 
{
    msg <- 'assay(<SummarizedExperiment>, i="numeric", ...) invalid subscript "i"'
    tryCatch({
        assays(x, ...)[[i]]
    }, error=function(err) {
        stop(msg, "\n", conditionMessage(err))
    })
})

setMethod(assay, c("SummarizedExperiment", "character"),
    function(x, i = names(x)[1], ...) 
{
    msg <- 'assay(<SummarizedExperiment>, i="character", ...) invalid subscript "i"'
    res <-
        tryCatch({
            assays(x, ...)[[i]]
        }, error=function(err) {
            stop(msg, "\n", conditionMessage(err))
        })
    if (is.null(res)) 
        stop(msg, "\n    '", i, "' not in names(assays(<SummarizedExperiment>))")
    res
})

setReplaceMethod("assay", c("SummarizedExperiment", "missing", "matrix"),
    function(x, i, ..., value)
{
    if (0L == length(assays(x)))
    {
        msg <- "'assay(<SummarizedExperiment>) <- value' length(assays(<SummarizedExperiment>)) is 0"
        stop(msg)
    }
    assays(x)[[1]] <- value
    x
})

setReplaceMethod("assay",
    c("SummarizedExperiment", "numeric", "matrix"),
    function(x, i = 1, ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

setReplaceMethod("assay",
    c("SummarizedExperiment", "character", "matrix"),
    function(x, i = names(x)[1], ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

## cannonical location for dim, dimnames; dimnames should be checked
## for consistency (if non-null) and stripped from assays on
## construction, or added from assays if row/col names are NULL in
## <SummarizedExperiment> but not assays. dimnames need to be added on
## to assays when assays() invoked
setMethod(dim, "SummarizedExperiment", function(x) {
        c(length(rowData(x)), nrow(colData(x)))
})

setMethod(dimnames, "SummarizedExperiment", function(x) {
        list(names(rowData(x)), rownames(colData(x)))
})

setReplaceMethod("dimnames", c("SummarizedExperiment", "list"),
    function(x, value) 
{
    rowData <- rowData(x)
    names(rowData) <- value[[1]]
    colData <- colData(x)
    rownames(colData) <- value[[2]]
    initialize(x, rowData=rowData, colData=colData)
})
 
setReplaceMethod("dimnames", c("SummarizedExperiment", "NULL"),
    function(x, value)
{
    callNextMethod(x, value=list(NULL, NULL))
})

## Subset -- array-like

.SummarizedExperiment.charbound <-
    function(idx, txt, fmt)
{
    orig <- idx
    idx <- match(idx, txt)
    if (any(bad <- is.na(idx))) {
        msg <- paste(IRanges:::selectSome(orig[bad]), collapse=" ")
        stop(sprintf(fmt, msg))
    }
    idx
}

.SummarizedExperiment.subset <-
    function(x, i, j, ...)
{
    if (is.character(i)) {
        msg <- "<SummarizedExperiment>[i,] index out of bounds: %s"
        i <- .SummarizedExperiment.charbound(i, rownames(x), msg)
    }
    if (is.character(j)) {
        msg <- "<SummarizedExperiment>[,j] index out of bounds: %s"
        j <- .SummarizedExperiment.charbound(j, colnames(x), msg)
    }

    assays <- endoapply(assays(x, withDimnames=FALSE), function(elt) {
                  if (class(elt) == "array")
                     elt[i, j, , drop=FALSE]
                  else
                     elt[i, j, drop=FALSE]
              })
    initialize(x, rowData=rowData(x)[i,,drop=FALSE],
               colData=colData(x)[j,,drop=FALSE],
               assays=assays)
}

setMethod("[", c("SummarizedExperiment", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,SummarizedExperiment,ANY,ANY-method'")
    if (missing(i) && missing(j))
        x
    else if (missing(i))
        .SummarizedExperiment.subset(x, TRUE, j, ...)
    else if (missing(j))
        .SummarizedExperiment.subset(x, i, TRUE, ...)
    else
        .SummarizedExperiment.subset(x, i, j, ...)
})

.SummarizedExperiment.subsetassign <-
    function(x, i, j, ..., value)
{
    if (is.character(i)) {
        msg <- "<SummarizedExperiment>[i,]<- index out of bounds: %s"
        i <- .SummarizedExperiment.charbound(i, rownames(x), msg)
    }
    if (is.character(j)) {
        msg <- "<SummarizedExperiment>[,j]<- index out of bounds: %s"
        j <- .SummarizedExperiment.charbound(j, colnames(x), msg)
    }
    initialize(x,
               exptData=c(exptData(x), exptData(value)),
               rowData=local({
                   r <- rowData(x)
                   r[i,] <- rowData(value)
                   names(r)[i] <- names(rowData(value))
                   r
               }), colData=local({
                   c <- colData(x)
                   c[j,] <- colData(value)
                   rownames(c)[j] <- rownames(colData(value))
                   c
               }), assays=local({
                   a <- assays(x, withDimnames=FALSE)
                   v <- assays(value, withDimnames=FALSE)
                   mendoapply(function(x, ..., value) {
                       x[i,j] <- value
                       x
                   }, x=a, value=v, ...)
               }))
}

setReplaceMethod("[",
    c("SummarizedExperiment", "ANY", "ANY", "SummarizedExperiment"),
    function(x, i, j, ..., value)
{
    if (missing(i) && missing(j))
        x
    else if (missing(i))
        .SummarizedExperiment.subsetassign(x, TRUE, j, ..., value=value)
    else if (missing(j))
        .SummarizedExperiment.subsetassign(x, i, TRUE, ..., value=value)
    else
        .SummarizedExperiment.subsetassign(x, i, j, ..., value=value)
})

## $, $<-, [[, [[<- for colData access

setMethod("[[", c("SummarizedExperiment", "ANY", "missing"),
    function(x, i, j, ...)
{
    colData(x)[[i, ...]]
})

setReplaceMethod("[[",
    c("SummarizedExperiment", "ANY", "missing", "ANY"),
    function(x, i, j, ..., value)
{
    colData(x)[[i, ...]] <- value
    x
})

setMethod("$", "SummarizedExperiment",
    function(x, name)
{
    colData(x)[[name]]
})

setReplaceMethod("$", c("SummarizedExperiment", "ANY"),
    function(x, name, value)
{
    colData(x)[[name]] <- value
    x
})

## show

setMethod(show, "SummarizedExperiment",
    function(object)
{
    selectSome <- IRanges:::selectSome
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(IRanges:::selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    expt <- names(exptData(object))
    if (is.null(expt))
        expt <- character(length(exptData(object)))
    scat("exptData(%d): %s\n", expt)
    nms <- names(assays(object, withDimnames=FALSE))
    if (is.null(nms))
        nms <- character(length(assays(object, withDimnames=FALSE)))
    scat("assays(%d): %s\n", nms)
    dimnames <- dimnames(object)
    dlen <- sapply(dimnames, length)
    if (dlen[[1]]) scat("rownames(%d): %s\n", dimnames[[1]])
    else scat("rownames: NULL\n")
    scat("rowData metadata column names(%d): %s\n",
         names(mcols(rowData(object))))
    if (dlen[[2]]) scat("colnames(%d): %s\n", dimnames[[2]])
    else cat("colnames: NULL\n")
    scat("colData names(%d): %s\n", names(colData(object)))
})
