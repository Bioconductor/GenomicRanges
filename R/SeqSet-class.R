setClass("SeqSet",
    representation(
        exptData="SimpleList",        # overall description
        rowData="GRanges",            # rows and their anntoations
        colData="DataFrame",          # columns and their annotations
        assays="SimpleList"))         # Data per-se -- e.g., list of matricies

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
        if (is.null(dim) || !all(xdim == edim))
            msg <- c(msg, sprintf(msg2, i))
    }

    msg
}

setValidity2("SeqSet", .valid.SeqSet)

## FIXME: MTM thinks that generics do not belong with class definitions

setGeneric("exptData", function(x, ...) standardGeneric("exptData"))

## rowData, colData seem too vague, but from eSet derived classes wanted to
## call the rows / cols something different from 'features' or 'samples', so
## might as well avoid the issue
setGeneric("rowData", function(x, ...) standardGeneric("rowData"))

setGeneric("colData", function(x, ...) standardGeneric("colData"))

setGeneric("assays", function(x, ...) standardGeneric("assays"))

setGeneric("assay", function(x, i, ...) standardGeneric("assay"))

## Simple 'getters'

setMethod(exptData, "SeqSet", function(x, ...) slot(x, "exptData"))

setMethod(rowData, "SeqSet", function(x, ...) slot(x, "rowData"))

setMethod(colData, "SeqSet", function(x, ...) slot(x, "colData"))

setMethod(assays, "SeqSet", function(x, ...) slot(x, "assays"))

## convenience for common use case 
setMethod(assay, c("SeqSet", "missing"),
    function(x, i, ...)
{
    assays <- assays(x)
    if (0L == length(assays))
    {
        msg <- 'assay(x="SeqSet", i="missing", ...) length(assays(x)) is 0'
        stop(msg)
    }
    assays[[1]]
})

setMethod(assay, c("SeqSet", "numeric"),
    function(x, i, ...) 
{
    msg <- 'assay(x="SeqSet", i="numeric", ...) invalid subscript "i"'
    tryCatch({
        assays(x)[[i]]
    }, error=function(err) {
        stop(msg, "\n", conditionMessage(err))
    })
})

setMethod(assay, c("SeqSet", "character"),
    function(x, i = names(x)[1], ...) 
{
    msg <- 'assay(x="SeqSet", i="character", ...) invalid subscript "i"'
    res <-
        tryCatch({
            assays(x)[[i]]
        }, error=function(err) {
            stop(msg, "\n", conditionMessage(err))
        })
    if (is.null(res))
        stop(msg, "\n    not in names(assays(x))")
    res
})

## cannonical location for dim, dimnames; dimnames should be checked for
## consistency (if non-null) and stripped from assays on construction, or
## added from assays if row/col names are NULL in <SeqSet> but not
## assays. dimnames need to be added on to assays whne assays() invoked
setMethod(dim, "SeqSet", function(x) {
    c(length(rowData(x)), length(colData(x)))
})

setMethod(dimnames, "SeqSet", function(x) {
    c(names(rowData(x)), names(colData(x)))
})
