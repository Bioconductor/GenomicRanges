setClass("SummarizedExperiment",
    representation(
        exptData="SimpleList",                # overall description
        rowData="GenomicRangesORGRangesList", # rows and their annotations
        colData="DataFrame",                  # columns and their annotations
        assays="SimpleList"),                 # Data -- e.g., list of matricies
    prototype(rowData=GRanges()))

.valid.SummarizedExperiment.assays_class <- function(x)
{
    ok <- sapply(assays(x, withDimnames=FALSE), function(cl) {
        is(cl, "matrix") || is(cl, "array") || is(cl, "Matrix")
    })
    if (!all(ok))
        return("'assays' must be class 'matrix', 'array', or 'Matrix'")
    NULL
}

.valid.SummarizedExperiment.rowData_dims <- function(x)
{
    if (!all(sapply(assays(x, withDimnames=FALSE), nrow) ==
             length(rowData(x))))
        return("'rowData' length differs from 'assays' nrow")
    NULL
}

.valid.SummarizedExperiment.colData_dims <- function(x)
{
    if (!all(sapply(assays(x, withDimnames=FALSE), ncol) ==
             nrow(colData(x))))
        return("'colData' nrow differs from 'assays' ncol")
    NULL
}

.valid.SummarizedExperiment.assays_dims <- function(x)
{
    c(.valid.SummarizedExperiment.rowData_dims(x),
      .valid.SummarizedExperiment.colData_dims(x))
}

.valid.SummarizedExperiment <- function(x)
{
    c(msg <- .valid.SummarizedExperiment.assays_class(x),
      if (is.null(msg)) {
          .valid.SummarizedExperiment.assays_dims(x)
      } else NULL)
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
            stop("'SummarizedExperiment' assay colnames must not be NULL")
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

## update / clone, see GenomicRanges-class.R

setMethod("clone", "SummarizedExperiment",
    function(x, ...)
{
    ## dput(slotNames("SummarizedExperiment"))
    valid_argnames <- c("exptData", "rowData", "colData", "assays")
    args <- IRanges:::extraArgsAsList(valid_argnames, ...)
    firstTime <- TRUE
    for (nm in names(args)) {
        if (firstTime) {
            slot(x, nm, FALSE) <- args[[nm]]
            firstTime <- FALSE
        }
        else {
            `slot<-`(x, nm, FALSE, args[[nm]])
        }
    }
    x
})

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
    clone(x, ..., exptData=value)
})

setReplaceMethod("exptData", c("SummarizedExperiment", "list"),
    function(x, ..., value)
{
    clone(x, ..., exptData=SimpleList(value))
})

setMethod(rowData, "SummarizedExperiment",
    function(x, ...) slot(x, "rowData"))

.SummarizedExperiment.rowData.replace <-
    function(x, ..., value)
{
    x <- clone(x, ..., rowData=value)
    msg <- .valid.SummarizedExperiment.rowData_dims(x)
    if (!is.null(msg))
        stop(msg)
    x
}

setReplaceMethod("rowData", c("SummarizedExperiment", "GenomicRanges"),
    .SummarizedExperiment.rowData.replace)

setReplaceMethod("rowData", c("SummarizedExperiment", "GRangesList"),
    .SummarizedExperiment.rowData.replace)

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
    clone(x, rowData=local({
        r <- rowData(x)
        mcols(r) <- value
        r
    }))
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
    x <- clone(x, ..., colData=value)
    msg <- .valid.SummarizedExperiment.colData_dims(x)
    if (!is.null(msg))
        stop(msg)
    x
})

setMethod(assays, "SummarizedExperiment",
    function(x, ..., withDimnames=TRUE)
{
    if (withDimnames)
        endoapply(slot(x, "assays"), "dimnames<-", dimnames(x))
    else
        slot(x, "assays")
})

.SummarizedExperiment.assays.replace <-
    function(x, ..., value)
{
    x <- clone(x, ..., assays=value)
    msg <- .valid.SummarizedExperiment(x)
    if (!is.null(msg))
        stop(msg)
    x
}

setReplaceMethod("assays", c("SummarizedExperiment", "SimpleList"),
    .SummarizedExperiment.assays.replace)

setReplaceMethod("assays", c("SummarizedExperiment", "list"),
    function(x, ..., value)
{
    value <- lapply(x, "dimnames<-", NULL)
    .SummarizedExperiment.assays.replace(x, ...,
                                         value=SimpleList(value))
})

## convenience for common use case 
setMethod(assay, c("SummarizedExperiment", "missing"),
    function(x, i, ...)
{
    assays <- assays(x, ...)
    if (0L == length(assays))
        stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
             "length(assays(<SummarizedExperiment>)) is 0'")
    assays[[1]]
})

setMethod(assay, c("SummarizedExperiment", "numeric"),
    function(x, i, ...) 
{
    tryCatch({
        assays(x, ...)[[i]]
    }, error=function(err) {
        stop("'assay(<", class(x), ">, i=\"numeric\", ...)' ",
             "invalid subscript 'i'\n", conditionMessage(err))
    })
})

setMethod(assay, c("SummarizedExperiment", "character"),
    function(x, i = names(x)[1], ...) 
{
    msg <- paste0("'assay(<", class(x), ">, i=\"character\", ...)' ",
                  "invalid subscript 'i'")
    res <- tryCatch({
        assays(x, ...)[[i]]
    }, error=function(err) {
        stop(msg, "\n", conditionMessage(err))
    })
    if (is.null(res))
        stop(msg, "\n'i' not in names(assays(<", class(x), ">))")
    res
})

setReplaceMethod("assay", c("SummarizedExperiment", "missing", "matrix"),
    function(x, i, ..., value)
{
    if (0L == length(assays(x)))
        stop("'assay(<", class(x), ">) <- value' ", "length(assays(<",
             class(x), ">)) is 0")
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
setMethod(dim, "SummarizedExperiment",
    function(x) 
{
    c(length(rowData(x)), nrow(colData(x)))
})

setMethod(dimnames, "SummarizedExperiment",
    function(x)
{
    list(names(rowData(x)), rownames(colData(x)))
})

setReplaceMethod("dimnames", c("SummarizedExperiment", "list"),
    function(x, value) 
{
    rowData <- rowData(x)
    names(rowData) <- value[[1]]
    colData <- colData(x)
    rownames(colData) <- value[[2]]
    clone(x, rowData=rowData, colData=colData)
})
 
setReplaceMethod("dimnames", c("SummarizedExperiment", "NULL"),
    function(x, value)
{
    dimnames(x) <- list(NULL, NULL)
    x
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

setMethod("[", c("SummarizedExperiment", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    if (missing(i) && missing(j))
        return(x)

    if (!missing(i) && is.character(i)) {
        fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
        i <- .SummarizedExperiment.charbound(i, rownames(x), fmt)
    }
    if (!missing(j) && is.character(j)) {
        fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
        j <- .SummarizedExperiment.charbound(j, colnames(x), fmt)
    }

    if (!missing(i) && !missing(j)) {
        fun <- function(elt)
        {
            if (length(dim(elt)) == 2L)
                elt[i, j, drop=FALSE]
            else
                elt[i, j, , drop=FALSE]
        }
        x <- clone(x, ..., rowData=rowData(x)[i],
            colData=colData(x)[j, , drop=FALSE],
            assays=endoapply(assays(x, withDimnames=FALSE), fun))
    } else if (missing(i)) {
        fun <- function(elt)
        {
            if (length(dim(elt)) == 2L)
                elt[, j, drop=FALSE]
            else
                elt[, j, , drop=FALSE]
        }
        x <- clone(x, ..., colData=colData(x)[j, , drop=FALSE],
            assays=endoapply(assays(x, withDimnames=FALSE), fun))
    } else {                            # missing(j)
        fun <- function(elt)
        {
            if (length(dim(elt)) == 2L)
                elt[i, , drop=FALSE]
            else
                elt[i, , , drop=FALSE]
        }
        x <- clone(x, ..., rowData=rowData(x)[i],
            assays=endoapply(assays(x, withDimnames=FALSE), fun))
    }
    x
})

setReplaceMethod("[",
    c("SummarizedExperiment", "ANY", "ANY", "SummarizedExperiment"),
    function(x, i, j, ..., value)
{
    if (missing(i) && missing(j))
        return(value)
    
    if (!missing(i) && is.character(i)) {
        fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
        i <- .SummarizedExperiment.charbound(i, rownames(x), fmt)
    }

    if (!missing(j) && is.character(j)) {
        fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
        j <- .SummarizedExperiment.charbound(j, colnames(x), fmt)
    }

    fun <- function(x, ..., value) {
        if (length(dim(x)) == 2)
            x[i, j] <- value
        else
            x[i, j, ] <- value
        x
    }

    if (!missing(i) && !missing(j)) {
        fun <- function(x, ..., value) {
            if (length(dim(x)) == 2)
                x[i, j] <- value
            else
                x[i, j, ] <- value
            x
        }
        x <- clone(x, ..., exptData=c(exptData(x), exptData(value)),
            rowData=local({
                r <- rowData(x)
                r[i] <- rowData(value)
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
                mendoapply(fun, x=a, value=v)
            }))
        msg <- .valid.SummarizedExperiment.assays_dims(x)
    } else if (missing(i)) {
        fun <- function(x, ..., value) {
            if (length(dim(x)) == 2)
                x[, j] <- value
            else
                x[, j, ] <- value
            x
        }
        x <- clone(x, ..., exptData=c(exptData(x), exptData(value)),
            colData=local({
                c <- colData(x)
                c[j,] <- colData(value)
                rownames(c)[j] <- rownames(colData(value))
                c
            }), assays=local({
                a <- assays(x, withDimnames=FALSE)
                v <- assays(value, withDimnames=FALSE)
                mendoapply(fun, x=a, value=v)
            }))
        msg <- .valid.SummarizedExperiment.colData_dims(x)
    } else {                            # missing(j)
        fun <- function(x, ..., value) {
            if (length(dim(x)) == 2)
                x[i, ] <- value
            else
                x[i, , ] <- value
            x
        }
        x <- clone(x, ..., exptData=c(exptData(x), exptData(value)),
            rowData=local({
                r <- rowData(x)
                r[i] <- rowData(value)
                names(r)[i] <- names(rowData(value))
                r
            }), assays=local({
                a <- assays(x, withDimnames=FALSE)
                v <- assays(value, withDimnames=FALSE)
                mendoapply(fun, x=a, value=v)
            }))
        msg <- .valid.SummarizedExperiment.rowData_dims(x)
    }
    if (!is.null(msg))
        stop(msg)
    x
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
