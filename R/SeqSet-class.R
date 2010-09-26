setClass("SeqSet",
    representation(
        exptData="SimpleList",        # overall description
        rowData="GenomicRanges",      # rows and their anntoations
        colData="DataFrame",          # columns and their annotations
        assays="SimpleList"),         # Data per-se -- e.g., list of matricies
    prototype(rowData=GRanges()))

.valid.SeqSet <- function(x)
{
    msg <- NULL
    msg1 <- '\n    is(assays(<SeqSet>)[[%d]], "matrix") is not TRUE'
    msg2 <- "\n    dim(<SeqSet>) does not equal dim(assays(<SeqSet>)[[%d]])"
    xdim <- dim(x)
    assays <- assays(x)
    for (i in seq_len(length(assays))) {
        if (!is(assays[[i]], "matrix")) {
            msg <- c(msg, sprintf(msg1, i))
            next
        }
        edim <- dim(assays[[i]])
        if (is.null(edim) || !all(xdim == edim))
            msg <- c(msg, sprintf(msg2, i))
    }

    msg
}

setValidity2("SeqSet", .valid.SeqSet)

SeqSet <-
    function(assays=SimpleList(), rowData=GRanges(),
             colData=DataFrame(), exptData=SimpleList(), ...,
             verbose=FALSE)
{
    ## FIXME: warn if dimnames(assays) != list( verbose=TRUE
    assays <- endoapply(assays, "dimnames<-", NULL)
    new("SeqSet", exptData=exptData, rowData=rowData, colData=colData,
        assays=assays, ...)
}

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

setGeneric("assays", function(x, ...) standardGeneric("assays"))

setGeneric("assays<-",
    function(x, ..., value) standardGeneric("assays<-"))

setGeneric("assay", function(x, i, ...) standardGeneric("assay"))

setGeneric("assay<-",
    function(x, i, ..., value) standardGeneric("assay<-"))

## Simple 'getters' / 'setters'

setMethod(exptData, "SeqSet", function(x, ...) slot(x, "exptData"))

setReplaceMethod("exptData", c("SeqSet", "SimpleList"),
    function(x, ..., value)
{
    initialize(x, ..., exptData=value)
})

setReplaceMethod("exptData", c("SeqSet", "list"),
    function(x, ..., value)
{
    initialize(x, ..., exptData=SimpleList(value))
})

setMethod(rowData, "SeqSet", function(x, ...) slot(x, "rowData"))

setReplaceMethod("rowData", c("SeqSet", "GenomicRanges"),
    function(x, ..., value)
{
    initialize(x, ..., rowData=value)
})

setMethod(colData, "SeqSet", function(x, ...) slot(x, "colData"))

setReplaceMethod("colData", c("SeqSet", "DataFrame"),
    function(x, ..., value)
{
    initialize(x, ..., colData=value)
})

setMethod(assays, "SeqSet",
    function(x, ...) 
{
    endoapply(slot(x, "assays"), "dimnames<-", dimnames(x))
})

setReplaceMethod("assays", c("SeqSet", "SimpleList"),
    function(x, ..., value)
{
    initialize(x, ..., assays=value)
})

setReplaceMethod("assays", c("SeqSet", "list"),
    function(x, ..., value)
{
    value <- lapply(x, "dimnames<-", NULL)
    initialize(x, ..., assays=SimpleList(value))
})

## convenience for common use case 
setMethod(assay, c("SeqSet", "missing"),
    function(x, i, ...)
{
    assays <- assays(x)
    if (0L == length(assays))
    {
        msg <- 'assay(<SeqSet>, i="missing", ...) length(assays(<SeqSet>)) is 0'
        stop(msg)
    }
    assays[[1]]
})

setMethod(assay, c("SeqSet", "numeric"),
    function(x, i, ...) 
{
    msg <- 'assay(<SeqSet>, i="numeric", ...) invalid subscript "i"'
    tryCatch({
        assays(x)[[i]]
    }, error=function(err) {
        stop(msg, "\n", conditionMessage(err))
    })
})

setMethod(assay, c("SeqSet", "character"),
    function(x, i = names(x)[1], ...) 
{
    msg <- 'assay(<SeqSet>, i="character", ...) invalid subscript "i"'
    res <-
        tryCatch({
            assays(x)[[i]]
        }, error=function(err) {
            stop(msg, "\n", conditionMessage(err))
        })
    if (is.null(res))
        stop(msg, "\n    '", i, "' not in names(assays(<SeqSet>))")
    res
})

setReplaceMethod("assay", c("SeqSet", "missing", "matrix"),
    function(x, i, ..., value)
{
    if (0L == length(assays(x)))
    {
        msg <- "'assay(<SeqSet>) <- value' length(assays(<SeqSet>)) is 0"
        stop(msg)
    }
    assays(x)[[1]] <- value
    x
})

setReplaceMethod("assay", c("SeqSet", "numeric", "matrix"),
    function(x, i = 1, ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

setReplaceMethod("assay", c("SeqSet", "character", "matrix"),
    function(x, i = names(x)[1], ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

## cannonical location for dim, dimnames; dimnames should be checked for
## consistency (if non-null) and stripped from assays on construction, or
## added from assays if row/col names are NULL in <SeqSet> but not
## assays. dimnames need to be added on to assays whne assays() invoked
setMethod(dim, "SeqSet", function(x) {
    c(length(rowData(x)), nrow(colData(x)))
})

setMethod(dimnames, "SeqSet", function(x) {
    list(names(rowData(x)), rownames(colData(x)))
})

setReplaceMethod("dimnames", c("SeqSet", "list"),
    function(x, value) 
{
    rowData <- rowData(x)
    names(rowData) <- value[[1]]
    colData <- colData(x)
    rownames(colData) <- value[[2]]
    initialize(x, rowData=rowData, colData=colData)
})
        
setReplaceMethod("dimnames", c("SeqSet", "NULL"),
    function(x, value)
{
    callNextMethod(x, value=list(NULL, NULL))
})

## Subset -- array-like; [[, $ not defined

.SeqSet.charbound <-
    function(idx, txt, fmt)
{
    orig <- idx
    idx <- match(idx, txt)
    if (any(bad <- is.na(idx))) {
        msg <- paste(selectSome(orig[bad]), collapse=" ")
        stop(sprintf(fmt, msg))
    }
    idx
}

.SeqSet.subset <-
    function(x, i, j, ...)
{
    if (is.character(i)) {
        msg <- "<SeqSet>[i,] index out of bounds: %s"
        i <- .SeqSet.charbound(i, rownames(x), msg)
    }
    if (is.character(j)) {
        msg <- "<SeqSet>[,j] index out of bounds: %s"
        j <- .SeqSet.charbound(j, colnames(x), msg)
    }
    initialize(x, rowData=rowData(x)[i,,drop=FALSE],
               colData=colData(x)[j,,drop=FALSE],
               assays=endoapply(assays(x), function(x) {
                   x <- x[i, j, drop=FALSE]
                   dimnames(x) <- NULL
                   x
               }))
}

setMethod("[", c("SeqSet", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,SeqSet,ANY,ANY-method'")
    if (missing(i) && missing(j))
        x
    else if (missing(i))
        .SeqSet.subset(x, TRUE, j, ...)
    else if (missing(j))
        .SeqSet.subset(x, i, TRUE, ...)
    else
        .SeqSet.subset(x, i, j, ...)
})

.SeqSet.subsetassign <-
    function(x, i, j, ..., value)
{
    if (is.character(i)) {
        msg <- "<SeqSet>[i,]<- index out of bounds: %s"
        i <- .SeqSet.charbound(i, rownames(x), msg)
    }
    if (is.character(j)) {
        msg <- "<SeqSet>[,j]<- index out of bounds: %s"
        j <- .SeqSet.charbound(j, colnames(x), msg)
    }
    initialize(x,
               exptData=c(exptData(x), exptData(value)),
               rowData=local({
                   r <- rowData(x)
                   r[i,] <- rowData(value)
                   r
               }), colData=local({
                   c <- colData(x)
                   c[j,] <- colData(value)
                   c
               }), assays=local({
                   a <- assays(x)
                   v <- assays(value)
                   mendoapply(function(x, ..., value) {
                       x[i,j] <- value
                       dimnames(x) <- NULL
                       x
                   }, x=a, value=v, ...)
               }))
}

setReplaceMethod("[", c("SeqSet", "ANY", "ANY", "SeqSet"),
    function(x, i, j, ..., value)
{
    if (missing(i) && missing(j))
        x
    else if (missing(i))
        .SeqSet.subsetassign(x, TRUE, j, ..., value=value)
    else if (missing(j))
        .SeqSet.subsetassign(x, i, TRUE, ..., value=value)
    else
        .SeqSet.subsetassign(x, i, j, ..., value=value)
})
           
setMethod(show, "SeqSet",
    function(object)
{
    selectSome <- IRanges:::selectSome
    scat <- function(fmt, vals)
    {
        ## FIXME: strwrap
        txt <- sprintf(fmt, length(vals),
                       paste(selectSome(vals), collapse=" "))
        cat(txt)
    }
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    scat("assays(%d): %s\n", names(assays(object)))
    dimnames <- dimnames(object)
    dlen <- sapply(dimnames, length)
    if (dlen[[1]]) scat("rownames(%d): %s\n", dimnames[[1]])
    else cat("rownames: NULL\n")
    if (dlen[[2]]) scat("colnames(%d): %s\n", dimnames[[2]])
    else cat("colnames: NULL\n")
})
