##
## Helper classes: Assays, ShallowData, ShallowSimpleListAssays

setClass("Assays")

.ShallowData <- setRefClass("ShallowData",
    fields = list( data = "ANY" ))

setMethod("clone", "ShallowData",  # not exported
    function(x, ...)
{
    args <- list(...)
    x <- x$copy()
    for (nm in names(args))
        x$field(nm, args[[nm]])
    x
})

.ShallowSimpleListAssays0 <- setRefClass("ShallowSimpleListAssays",
    fields = list( data = "SimpleList" ),
    contains = c("ShallowData", "Assays"))

.ShallowSimpleListAssays <-
    function(..., data = SimpleList())
    ## this two-step constructor is much faster, 16 May, 2013
{
    xx <- .ShallowSimpleListAssays0(...)
    xx$data <- data
    xx
}

##
## SummarizedExperiment

setClass("SummarizedExperiment",
    representation(
        exptData="SimpleList",                # overall description
        rowData="GenomicRangesORGRangesList", # rows and their annotations
        colData="DataFrame",                  # columns and their annotations
        assays="Assays"),                     # Data -- e.g., list of matricies
    prototype(
        rowData=GRanges(),
        assays=.ShallowSimpleListAssays(data=SimpleList())))

## validity

.valid.SummarizedExperiment.assays_current <- function(x)
{
    if (!is(slot(x, "assays"), "Assays"))
        return("'assays' is out-of-date; use updateObject()")
    NULL
}

.valid.SummarizedExperiment.assays_class <- function(x)
{
    ok <- sapply(assays(x, withDimnames=FALSE), function(cl) {
        (!is.null(dim(cl))) && (length(dim(cl)) >= 2L)
    })
    if (!all(ok))
        return("'assays' must be matrix-like with 2 (or more?) dimensions")
    NULL
}

.valid.SummarizedExperiment.rowRanges_dims <- function(x)
{
    if (!all(sapply(assays(x, withDimnames=FALSE), nrow) ==
             length(rowRanges(x))))
        return("'rowRanges' length differs from 'assays' nrow")
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
    c(.valid.SummarizedExperiment.rowRanges_dims(x),
      .valid.SummarizedExperiment.colData_dims(x))
}

.valid.SummarizedExperiment <- function(x)
{
    c(.valid.SummarizedExperiment.assays_current(x),
      msg <- .valid.SummarizedExperiment.assays_class(x),
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
   function(assays, rowRanges=GRangesList(), colData=DataFrame(),
             exptData=SimpleList(), ...,
             verbose=FALSE)
{
    extra_args <- list(...)
    rowData <- extra_args$rowData
    if (!is.null(rowData)) {
        if (!missing(rowRanges))
            stop("both 'rowRanges' and 'rowData' are specified")
        msg <- c("The 'rowData' argument is deprecated. ",
                 "Please use 'rowRanges' instead.")
        .Deprecated(msg=msg)
        rowRanges <- rowData
        extra_args$rowData <- NULL
    } else if (missing(rowRanges) && 0L != length(assays)) {
        rowRanges <- .GRangesList_assays(assays)
    }
    if (missing(colData) && 0L != length(assays)) {
        nms <- colnames(assays[[1]])
        if (is.null(nms) && 0L != ncol(assays[[1]]))
            stop("'SummarizedExperiment' assay colnames must not be NULL")
        colData <- DataFrame(row.names=nms)
    }

    FUN <- function(x) {
        exp <- list(names(rowRanges), rownames(colData))
        ## dimnames as NULL or list(NULL, NULL)
        all(sapply(dimnames(x), is.null)) ||
            ## or consistent with row / colData
            identical(dimnames(x)[1:2], exp)
    }
    if (!all(sapply(assays, FUN)))
        assays <- endoapply(assays, unname)
    if (!is(assays, "Assays"))
        assays <- .ShallowSimpleListAssays(data=assays)
    do.call(new, c(list("SummarizedExperiment",
                        exptData=exptData,
                        rowData=rowRanges,
                        colData=colData,
                        assays=assays),
                   extra_args))
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
    if (is.list(assays))
        ## special case -- matrix of lists
        assays <- list(assays)
    SummarizedExperiment(SimpleList(assays), ...)
})

## update / clone

setMethod("clone", "SummarizedExperiment",  # not exported
    function(x, ...)
{
    ## S4Vectors:::extraArgsAsList would prevent using clone on
    ## subclasses
    args <- list(...)
    firstTime <- TRUE
    for (nm in names(args)) {
        s <- slot(x, nm)
        v <- args[[nm]]
        if (is(s, "ShallowData"))
            v <- clone(s, data=v)
        if (firstTime) {
            slot(x, nm, FALSE) <- v
            firstTime <- FALSE
        } else {
            `slot<-`(x, nm, FALSE, v)
        }
    }
    x
})

setGeneric("value",  # not exported
    function(x, name, ...) standardGeneric("value"),
    signature = "x")

setMethod("value", "SummarizedExperiment", # not exported
    function(x, name, ...)
{
    s <- slot(x, name)
    if (is(s, "ShallowData"))
        s <- s$data
    s
})

## getter / setter generics

setGeneric("exptData", function(x, ...) standardGeneric("exptData"))

setGeneric("exptData<-",
    function(x, ..., value) standardGeneric("exptData<-"))

## rowRanges, colData seem too vague, but from eSet derived classes wanted to
## call the rows / cols something different from 'features' or 'samples', so
## might as well avoid the issue
setGeneric("rowRanges", function(x, ...) standardGeneric("rowRanges"))

setGeneric("rowRanges<-",
    function(x, ..., value) standardGeneric("rowRanges<-"))

rowData <- function(...)
{
    .Deprecated("rowRanges")
    rowRanges(...)
}

`rowData<-` <- function(x, ..., value)
{
    .Deprecated("rowRanges<-")
    rowRanges(x, ...) <- value
}

setGeneric("colData", function(x, ...) standardGeneric("colData"))

setGeneric("colData<-",
    function(x, ..., value) standardGeneric("colData<-"))

setGeneric("assayNames", function(x, ...) standardGeneric("assayNames"))

setGeneric("assayNames<-",
    function(x, ..., value) standardGeneric("assayNames<-"))

setGeneric("assays",
    function(x, ..., withDimnames=TRUE) standardGeneric("assays"),
    signature="x")

setGeneric("assays<-",
    function(x, ..., withDimnames=TRUE, value) standardGeneric("assays<-"),
    signature=c("x", "value"))

setGeneric("assay", function(x, i, ...) standardGeneric("assay"))

setGeneric("assay<-",
    function(x, i, ..., value) standardGeneric("assay<-"))

## Simple 'getters' / 'setters'

setMethod(exptData, "SummarizedExperiment",
    function(x, ...) value(x, "exptData"))

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

setMethod(rowRanges, "SummarizedExperiment",
    function(x, ...) value(x, "rowData"))

.SummarizedExperiment.rowRanges.replace <-
    function(x, ..., value)
{
    x <- clone(x, ..., rowData=value)
    msg <- .valid.SummarizedExperiment.rowRanges_dims(x)
    if (!is.null(msg))
        stop(msg)
    x
}

setReplaceMethod("rowRanges", c("SummarizedExperiment", "GenomicRanges"),
    .SummarizedExperiment.rowRanges.replace)

setReplaceMethod("rowRanges", c("SummarizedExperiment", "GRangesList"),
    .SummarizedExperiment.rowRanges.replace)

setMethod(colData, "SummarizedExperiment",
    function(x, ...) value(x, "colData"))

setReplaceMethod("colData", c("SummarizedExperiment", "DataFrame"),
    function(x, ..., value)
{
    x <- clone(x, ..., colData=value)
    msg <- .valid.SummarizedExperiment.colData_dims(x)
    if (!is.null(msg))
        stop(msg)
    x
})

setMethod("assayNames", "SummarizedExperiment",
    function(x, ...)
{
    names(assays(x, withDimnames=FALSE))
})

setMethod("assayNames<-", c("SummarizedExperiment", "character"),
    function(x, ..., value)
{
    names(assays(x, withDimnames=FALSE)) <- value
    x
})

setMethod(assays, "SummarizedExperiment",
    function(x, ..., withDimnames=TRUE)
{
    if (withDimnames)
        endoapply(value(x, "assays"), "dimnames<-", dimnames(x))
    else
        value(x, "assays")
})

.SummarizedExperiment.assays.replace <-
    function(x, ..., withDimnames=TRUE, value)
{
    ## withDimnames arg allows names(assays(se, withDimnames=FALSE)) <- value
    ok <- vapply(value, function(elt, xdimnames) {
        e <- dimnames(elt)
        (is.null(e[[1]]) || identical(e[[1]], xdimnames[[1]])) &&
             (is.null(e[[2]]) || identical(e[[2]], xdimnames[[2]]))
    }, logical(1), xdimnames=dimnames(x))
    if (!all(ok))
        stop("current and replacement dimnames() differ")
    x <- clone(x, ..., assays=value)
    msg <- .valid.SummarizedExperiment(x)
    if (!is.null(msg))
        stop(msg)
    x
}

setReplaceMethod("assays", c("SummarizedExperiment", "SimpleList"),
    .SummarizedExperiment.assays.replace)

setReplaceMethod("assays", c("SummarizedExperiment", "list"),
    .SummarizedExperiment.assays.replace)

## convenience for common use case
setMethod(assay, c("SummarizedExperiment", "missing"),
    function(x, i, ...)
{
    assays <- assays(x, ...)
    if (0L == length(assays))
        stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
             "length(assays(<", class(x), ">)) is 0'")
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
    c(length(rowRanges(x)), nrow(colData(x)))
})

setMethod(dimnames, "SummarizedExperiment",
    function(x)
{
    list(names(rowRanges(x)), rownames(colData(x)))
})

setReplaceMethod("dimnames", c("SummarizedExperiment", "list"),
    function(x, value)
{
    rowRanges <- rowRanges(x)
    names(rowRanges) <- value[[1]]
    colData <- colData(x)
    rownames(colData) <- value[[2]]
    clone(x, rowData=rowRanges, colData=colData)
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
        msg <- paste(BiocGenerics:::selectSome(orig[bad]), collapse=" ")
        stop(sprintf(fmt, msg))
    }
    idx
}

.SummarizedExperiment.assays.subset <-
    function(x, i, j, ...)
{
    ## need to expand Rle's for subsetting standard matrix
    if (!missing(i) && !missing(j)) {
        fun <- function(x) {
            switch(length(dim(x)),
                   stop("'[' on assays() with 1 dimension not supported"),
                   x[i, j, drop=FALSE],
                   x[i, j, , drop=FALSE],
                   x[i, j, , , drop=FALSE],
                   stop("'[' on assays() with >4 dimensions not supported"))
        }
    } else if (!missing(i)) {
        fun <- function(x) {
            switch(length(dim(x)),
                   stop("'[' on assays() with 1 dimension not supported"),
                   x[i, , drop=FALSE],
                   x[i, , , drop=FALSE],
                   x[i, , , , drop=FALSE],
                   stop("'[' on assays() with >4 dimensions not supported"))
        }
    } else if (!missing(j)) {
        fun <- function(x) {
            switch(length(dim(x)),
                   stop("'[' on assays() with 1 dimension not supported"),
                   x[, j, drop=FALSE],
                   x[, j, , drop=FALSE],
                   x[, j, , , drop=FALSE],
                   stop("'[' on assays() with >4 dimensions not supported"))
        }
    }
    endoapply(assays(x, withDimnames=FALSE), fun)
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
        ii <- as.vector(i)
        jj <- as.vector(j)
        x <- clone(x, ..., rowData=rowRanges(x)[i],
            colData=colData(x)[j, , drop=FALSE],
            assays=.SummarizedExperiment.assays.subset(x, ii, jj))
    } else if (missing(i)) {
        jj <- as.vector(j)
        x <- clone(x, ..., colData=colData(x)[j, , drop=FALSE],
            assays=.SummarizedExperiment.assays.subset(x, j=jj))
    } else {                            # missing(j)
        ii <- as.vector(i)
        x <- clone(x, ..., rowData=rowRanges(x)[i],
            assays=.SummarizedExperiment.assays.subset(x, ii))
    }
    x
})

.SummarizedExperiment.assays.subsetgets <-
    function(x, i, j, ..., value)
{
    ## need to expand Rle's for subsetting standard matrix
    if (!missing(i) && !missing(j)) {
        fun <- function(x, value) {
            switch(length(dim(x)),
                   stop("'[<-' on assays() with 1 dimension not supported"),
                   x[i, j] <- value,
                   x[i, j, ] <- value,
                   x[i, j, , ] <- value,
                   stop("'[<-' on assays() with >4 dimensions not supported"))
            x
        }
    } else if (!missing(i)) {
        fun <- function(x, value) {
            switch(length(dim(x)),
                   stop("'[<-' on assays() with 1 dimension not supported"),
                   x[i, ] <- value,
                   x[i, , ] <- value,
                   x[i, , , ] <- value,
                   stop("'[<-' on assays() with >4 dimensions not supported"))
            x
        }
    } else if (!missing(j)) {
        fun <- function(x, value) {
            switch(length(dim(x)),
                   stop("'[<-' on assays() with 1 dimension not supported"),
                   x[, j] <- value,
                   x[, j, ] <- value,
                   x[, j, , ] <- value,
                   stop("'[<-' on assays() with >4 dimensions not supported"))
            x
        }
    }
    a <- assays(x, withDimnames=FALSE)
    v <- assays(value, withDimnames=FALSE)
    mendoapply(fun, x=a, value=v)
}

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

    if (!missing(i) && !missing(j)) {
        ii <- as.vector(i)
        jj <- as.vector(j)
        x <- clone(x, ..., exptData=c(exptData(x), exptData(value)),
            rowData=local({
                r <- rowRanges(x)
                r[i] <- rowRanges(value)
                names(r)[ii] <- names(rowRanges(value))
                r
            }), colData=local({
                c <- colData(x)
                c[j,] <- colData(value)
                rownames(c)[jj] <- rownames(colData(value))
                c
            }), assays=.SummarizedExperiment.assays.subsetgets(x, ii, jj,
                  ..., value=value))
        msg <- .valid.SummarizedExperiment.assays_dims(x)
    } else if (missing(i)) {
        jj <- as.vector(j)
        x <- clone(x, ..., exptData=c(exptData(x), exptData(value)),
            colData=local({
                c <- colData(x)
                c[j,] <- colData(value)
                rownames(c)[jj] <- rownames(colData(value))
                c
            }), assays=.SummarizedExperiment.assays.subsetgets(x, j=jj,
                  ..., value=value))
        msg <- .valid.SummarizedExperiment.colData_dims(x)
    } else {                            # missing(j)
        ii <- as.vector(i)
        x <- clone(x, ..., exptData=c(exptData(x), exptData(value)),
            rowData=local({
                r <- rowRanges(x)
                r[i] <- rowRanges(value)
                names(r)[ii] <- names(rowRanges(value))
                r
            }), assays=.SummarizedExperiment.assays.subsetgets(x, ii,
                  ..., value=value))
        msg <- .valid.SummarizedExperiment.rowRanges_dims(x)
    }
    if (!is.null(msg))
        stop(msg)
    x
})

## rbind, cbind

## Appropriate for objects with different ranges and same samples.
setMethod("rbind", "SummarizedExperiment",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    .rbind.SummarizedExperiment(args)
})

.rbind.SummarizedExperiment <- function(args)
{
    if (!.compare(lapply(args, colnames)))
            stop("'...' objects must have the same colnames")
    if (!.compare(lapply(args, ncol)))
            stop("'...' objects must have the same number of samples")
    rowRanges <- do.call(c, lapply(args, rowRanges))
    colData <- .cbind.DataFrame(args, colData, "colData")
    assays <- .bind.arrays(args, rbind, "assays")
    exptData <- do.call(c, lapply(args, exptData))

    initialize(args[[1]], assays=assays, rowData=rowRanges,
               colData=colData, exptData=exptData)
}

## Appropriate for objects with same ranges and different samples.
setMethod("cbind", "SummarizedExperiment",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    .cbind.SummarizedExperiment(args)
})

.cbind.SummarizedExperiment <- function(args)
{
    if (!.compare(lapply(args, rowRanges), TRUE))
        stop("'...' object ranges (rows) are not compatible")
    rowRanges <- rowRanges(args[[1]])
    mcols(rowRanges) <- .cbind.DataFrame(args, mcols, "mcols")
    colData <- do.call(rbind, lapply(args, colData))
    assays <- .bind.arrays(args, cbind, "assays")
    exptData <- do.call(c, lapply(args, exptData))

    initialize(args[[1]], assays=assays, rowData=rowRanges,
               colData=colData, exptData=exptData)
}

.compare <- function(x, GenomicRanges=FALSE)
{
    x1 <- x[[1]]
    if (GenomicRanges) {
        if (is(x1, "GRangesList")) {
            x <- lapply(x, unlist)
            x1 <- x[[1]]
        }
        for (i in seq_along(x)[-1]) {
            if (length(x1) != length(x[[i]]))
                return(FALSE)
            ok <- x1 == x[[i]]
            if (!all(ok))
                return(FALSE)
        }
        return(TRUE)
    } else {
        all(sapply(x[-1],
            function(xelt) all(identical(xelt, x[[1]]))))
    }
}

.cbind.DataFrame <- function(args, accessor, accessorName)
{
    lst <- lapply(args, accessor)
    if (!.compare(lst)) {
        nms <- lapply(lst, names)
        nmsv <- unlist(nms, use.names=FALSE)
        names(nmsv) <- rep(seq_along(nms), elementLengths(nms))
        dups <- duplicated(nmsv)
        ## no duplicates
        if (!any(dups))
            return(do.call(cbind, lst))
        ## confirm duplicates are the same
        lapply(nmsv[duplicated(nmsv)], function(d) {
            if (!.compare(lapply(lst, "[", d)))
                stop("column(s) '", unname(d),
                     "' in ", sQuote(accessorName),
                     " are duplicated and the data do not match")})
        ## remove duplicates
        do.call(cbind, lst)[,!dups]
    } else {
        lst[[1]]
    }
}

.bind.array.elements <- function(index, lst, bind) {
    e1 <- lapply(lst, "[[", index)
    dim <- .get.assay.dimension(e1, bind)
    if (is.na(dim[3])) {
        do.call(bind, e1)
    } else {
        e2 <- lapply(1:dim[3], function(i) {
            do.call(bind, lapply(e1, "[", ,,i))
        })
        array(do.call(c, e2), dim=dim)
    }
}

.bind.arrays <- function(args, bind, accessor)
{
    lst <- lapply(args, accessor)
    if (!length(elts <- unique(elementLengths(lst))) == 1L)
        stop("elements in ", sQuote(accessor),
             " must have the same length")
    if (elts == 0L)
        return(.ShallowSimpleListAssays(data=SimpleList()))
    var <- lapply(lst,  names)
    if (is.null(uvar <- unique(unlist(var)))) {
        ## no names, match by position
        res <- lapply(seq_along(elts), .bind.array.elements, lst=lst, bind=bind)
    } else {
        ## match by name
        if (!.compare(var))
            stop("elements in ", sQuote(accessor),
                 " must have the same names")
        res <- lapply(uvar, .bind.array.elements, lst=lst, bind=bind)
        names(res) <- uvar
    }
    .ShallowSimpleListAssays(data=SimpleList(res))
}

.get.assay.dimension <- function(lst, bind)
{
    dim <- lapply(lst, dim)
    if (!.compare(lapply(dim, "[", 3)))
        stop("elements in assays must have the same dimension")
    if (identical(bind, cbind))
        c(dim[[1]][1], do.call(sum, lapply(dim, "[", 2)), dim[[1]][3])
    else
        c(do.call(sum, lapply(dim, "[", 1)), dim[[1]][2], dim[[1]][3])
}

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

.DollarNames.SummarizedExperiment <- function(x, pattern)
    grep(pattern, names(colData(x)), value=TRUE)

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

## compatibility

setMethod(updateObject, "SummarizedExperiment",
    function(object, ..., verbose=FALSE)
{
    s <- slot(object, "assays")
    if (is(s, "SimpleList"))
        slot(object, "assays") <- .ShallowSimpleListAssays(data=s)
    object
})

## show

setMethod(show, "SummarizedExperiment",
    function(object)
{
    selectSome <- BiocGenerics:::selectSome
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(BiocGenerics:::selectSome(vals), collapse=" ")
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
    scat("rowRanges metadata column names(%d): %s\n",
         names(mcols(rowRanges(object))))
    if (dlen[[2]]) scat("colnames(%d): %s\n", dimnames[[2]])
    else cat("colnames: NULL\n")
    scat("colData names(%d): %s\n", names(colData(object)))
})

.from_rowRanges_to_FeatureData <- function(from)
{
    if (is(from, "GRanges")) {
        fd <- .from_GRanges_to_FeatureData(from)
    } else if (is(from, "GRangesList")) {
        fd <- .from_GRangesList_to_FeatureData(from)
    } else {
        stop("class: ", sQuote(class(from)), " is not a supported type for rowRanges coercion")
    }
    Biobase::featureNames(fd) <- names(from)
    fd
}

.from_GRanges_to_FeatureData <- function(from)
{
    data <- as.data.frame(from)

    ## the first mcols are automatically included in the data.frame from
    ## as.data.frame, the secondary mcols holds the metadata for the first
    ## metadata columns.
    metaData <- mcols(mcols(from))
    if (is.null(metaData)) {
        metaData <- as.data.frame(matrix(ncol=0, nrow=NCOL(data)))
    } else {
        metaData <- as.data.frame(metaData)
    }
    Biobase::AnnotatedDataFrame(data, metaData)
}
.from_GRangesList_to_FeatureData <- function(from)
{
    data <- as.data.frame(mcols(from))

    ## the first mcols are automatically included in the data.frame from
    ## as.data.frame, the secondary mcols holds the metadata for the first
    ## metadata columns.
    metaData <- mcols(mcols(from))
    if (is.null(metaData)) {
        metaData <- as.data.frame(matrix(ncol=0, nrow=NCOL(data)))
    } else {
        metaData <- as.data.frame(metaData)
    }
    Biobase::AnnotatedDataFrame(data, metaData)
}

.from_AnnotatedDataFrame_to_DataFrame <- function(from)
{
    df <- DataFrame(Biobase::pData(from))
    mcols(df) <- DataFrame(Biobase::varMetadata(from))
    df
}

.from_DataFrame_to_AnnotatedDataFrame <- function(df)
{
    data <- as(df, "data.frame")
    metaData <- mcols(df)
    if (is.null(metaData)) {
        metaData <- as.data.frame(matrix(ncol=0, nrow=NCOL(data)))
    } else {
        metaData <- as(metaData, "data.frame")
    }
    Biobase::AnnotatedDataFrame(data, metaData)
}

# If the ExpressionSet has featureData with range information make GRanges out
# of that, otherwise make an empty GRangesList with names from the featureNames
naiveRangeMapper <- function(from)
{
    nms <- Biobase::featureNames(from)
    res <- tryCatch(makeGRangesFromDataFrame(Biobase::pData(Biobase::featureData(from)), keep.extra.columns = TRUE),
             error = function(e) {
                 res <- relist(GRanges(), vector("list", length=length(nms)))
                 mcols(res) <- .from_AnnotatedDataFrame_to_DataFrame(Biobase::featureData(from))
                 res
             })
    names(res) <- nms
    res
}

# Simple ProbeId to Range mapper
# Probes with multiple ranges are dropped
# The sign of the chromosome location is assumed to contain the strand
# information
probeRangeMapper <- function(from)
{
    annotation <- Biobase::annotation(from)
    if (identical(annotation, character(0))) {
        return(naiveRangeMapper(from))
    }
    if (requireNamespace("annotate", quietly = TRUE)) {
        annotationPackage <- annotate::annPkgName(annotation)
        if (require(annotationPackage, character.only = TRUE, quietly = TRUE)) {
            db <- get(annotationPackage, envir = asNamespace(annotationPackage))
            pid <- Biobase::featureNames(from)
            locs <- AnnotationDbi::select(db, pid, columns = c("CHR", "CHRLOC", "CHRLOCEND"))
            locs <- na.omit(locs)
            dups <- duplicated(locs$PROBEID)
            if (any(dups)) {
                locs <- locs[!dups, , drop = FALSE]
            }
            strand <- ifelse(locs$CHRLOC > 0, "+", "-")
            res <- GRanges(seqnames = locs$CHR,
                           ranges = IRanges(abs(locs$CHRLOC), abs(locs$CHRLOCEND)),
                           strand = strand)
            names(res) <- locs$PROBEID

            if (NROW(res) < length(pid)) {
                warning(length(pid) - NROW(res), " probes could not be mapped.", call. = FALSE)
            }
            res
        } else {
            stop("Failed to load ", sQuote(annotationPackage), " package", call. = FALSE)
        }
    } else {
        stop("Failed to load annotate package", call. = FALSE)
    }
}

# Simple ProbeId to Gene mapper
# Is there a way to get the txDb given the annotation package?
geneRangeMapper <- function(txDbPackage, key = "ENTREZID")
{
    function(from) {
        annotation <- Biobase::annotation(from)
        if (identical(annotation, character(0))) {
            return(naiveRangeMapper(from))
        }
        if (requireNamespace("annotate", quietly = TRUE)) {
            annotationPackage <- annotate::annPkgName(annotation)
            if (require(annotationPackage, character.only = TRUE, quietly = TRUE)) {
                db <- get(annotationPackage, envir = asNamespace(annotationPackage))
                pid <- Biobase::featureNames(from)
                probeIdToGeneId <- AnnotationDbi::mapIds(db, pid, key, "PROBEID")
                geneIdToProbeId <- setNames(names(probeIdToGeneId), probeIdToGeneId)

                if (requireNamespace(txDbPackage, quietly = TRUE)) {
                    txDb <- get(txDbPackage, envir = asNamespace(txDbPackage))
                    genes <- GenomicFeatures::genes(txDb)
                    probesWithAMatch <- probeIdToGeneId[probeIdToGeneId %in% names(genes)]
                    res <- genes[probesWithAMatch]
                    names(res) <- geneIdToProbeId[names(res)]
                    if (NROW(res) < length(pid)) {
                        warning(length(pid) - NROW(res), " probes could not be mapped.", call. = FALSE)
                    }
                    res
                } else {
                    stop("Failed to load ", sQuote(txDbPackage), " package", call. = FALSE)
                }

            } else {
                stop("Failed to load ", sQuote(annotationPackage), " package", call. = FALSE)
            }
        } else {
            stop("Failed to load annotate package", call. = FALSE)
        }
    }
}

makeSummarizedExperimentFromExpressionSet <- function(from, mapFun = naiveRangeMapper)
{
    if (requireNamespace("Biobase", quietly = TRUE)) {
        mapFun <- match.fun(mapFun)
        rowRanges <- mapFun(from)
        matches <- match(names(rowRanges),
                         Biobase::featureNames(from),
                         nomatch = 0)
        from <- from[matches, drop = FALSE]
        assays <- as.list(Biobase::assayData(from))
        colData <- .from_AnnotatedDataFrame_to_DataFrame(Biobase::phenoData(from))
        exptData <- SimpleList(
                               experimentData = Biobase::experimentData(from),
                               annotation = Biobase::annotation(from),
                               protocolData = Biobase::protocolData(from)
                               )

        SummarizedExperiment(
                             assays = assays,
                             rowRanges = rowRanges,
                             colData = colData,
                             exptData = exptData
                             )
    }
}

suppressMessages(
    setAs("ExpressionSet", "SummarizedExperiment", function(from)
    {
        makeSummarizedExperimentFromExpressionSet(from, naiveRangeMapper)
    })
)

suppressMessages(
    setAs("SummarizedExperiment", "ExpressionSet", function(from)
    {
        if (requireNamespace("Biobase", quietly = TRUE)) {

            assayData <- list2env(as.list(assays(from)))

            numAssays <- length(assayData)

            if(numAssays == 0) {
                assayData$exprs <- new("matrix")
            } else if (!"exprs" %in% ls(assayData)) {
                ## if there isn't an exprs assay we need to pick one as exprs, so
                ## rename the first element exprs and issue a warning.
                exprs <- ls(assayData)[[1]]
                warning("No assay named ", sQuote("exprs"), " found, renaming ",
                        exprs, " to ", sQuote("exprs"), ".")
                assayData[["exprs"]] <- assayData[[exprs]]
                rm(list=exprs, envir=assayData)
            }

            featureData <- .from_rowRanges_to_FeatureData(rowRanges(from))
            phenoData <- .from_DataFrame_to_AnnotatedDataFrame(colData(from))

            ed <- exptData(from)

            experimentData <- if (!is.null(ed$experimentData)) {
                ed$experimentData
            } else {
                Biobase::MIAME()
            }

            annotation <- if (!is.null(ed$annotation)) {
                ed$annotation
            } else {
                character()
            }

            protocolData <- if (!is.null(ed$protocolData)) {
                ed$protocolData
            } else {
                Biobase::annotatedDataFrameFrom(assayData, byrow=FALSE)
            }

            Biobase::ExpressionSet(assayData,
                                   phenoData = phenoData,
                                   featureData = featureData,
                                   experimentData = experimentData,
                                   annotation = annotation,
                                   protocolData = protocolData
                                   )
        }
    })
)
