## colData-as-GRanges compatibility: allow direct access to GRanges /
## GRangesList colData for select functions

## Not supported:
## 
## Not consistent SummarizedExperiment structure: length, names,
##   as.data.frame, c, splitAsListReturnedClass.
## Length-changing endomorphisms: disjoin, gaps, reduce,
##   unique.
## 'legacy' data types / functions: as "RangedData", as "RangesList",
##   renameSeqlevels, keepSeqlevels.
## Possile to implement, but not yet: Ops, resolveHits, map,
##   seqselect, seqselect<-, window, window<-

## mcols
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

## Single dispatch, generic signature fun(x, ...)
local({
    .funs <-
        c("coverage", "disjointBins", "duplicated", "end", "end<-",
          "flank", "isDisjoint", "narrow", "ranges", "resize",
          "restrict", "seqinfo", "seqnames", "shift", "start",
          "start<-", "strand", "width", "width<-")

    endomorphisms <-
        c(.funs[grepl("<-$", .funs)], "flank", "narrow", "resize",
          "restrict", "shift")

    tmpl <- function() {}
    environment(tmpl) <- parent.frame(2)
    for (.fun in .funs) {
        generic <- getGeneric(.fun)
        formals(tmpl) <- formals(generic)
        fmls <- as.list(formals(tmpl))
        fmls[] <- sapply(names(fmls), as.symbol)
        fmls[[generic@signature]] <- quote(rowData(x))
        if (.fun %in% endomorphisms)
            body(tmpl) <- substitute({
                rowData(x) <- do.call(FUN, ARGS)
                x
            }, list(FUN=.fun, ARGS=fmls))
        else
            body(tmpl) <-
                substitute(do.call(FUN, ARGS),
                           list(FUN=as.symbol(.fun), ARGS=fmls))
        setMethod(.fun, "SummarizedExperiment", tmpl)
    }
})

setMethod("granges", "SummarizedExperiment",
    function(x, ...) rowData(x))

## 2-argument dispatch:
## compare / Compare 
## precede, follow, nearest, distance, distanceToNearest
## 
## findOverlaps / countOverlaps, match, subsetByOverlaps: see
## findOverlaps-method.R
.SummarizedExperiment.compare <-
    function(x, y)
{
    if (is(x, "SummarizedExperiment"))
        x <- rowData(x)
    if (is(y, "SummarizedExperiment"))
        y <- rowData(y)
    compare(x, y)
}

.SummarizedExperiment.Compare <-
    function(e1, e2)
{
    if (is(e1, "SummarizedExperiment"))
        e1 <- rowData(e1)
    if (is(e2, "SummarizedExperiment"))
        e2 <- rowData(e2)
    callGeneric(e1=e1, e2=e2)
}

.SummarizedExperiment.nearest.missing <-
    function(x, subject, select = c("arbitrary", "all"),
             ignore.strand = FALSE, ...)
{
    select <- match.arg(select)
    x <- rowData(x)
    nearest(x=x, subject=x, select=select,
            ignore.strand=ignore.strand, ignoreSelf=TRUE, ...)
}

.SummarizedExperiment.distance <-
    function(x, y, ignore.strand = FALSE, ...)
{
    if (is(x, "SummarizedExperiment"))
        x <- rowData(x)
    if (is(y, "SummarizedExperiment"))
        y <- rowData(y)
    distance(x, y, ignore.strand=ignore.strand, ...)
}

.SummarizedExperiment.distanceToNearest <-
    function(x, subject, ignore.strand = FALSE, ...)
{
    if (is(x, "SummarizedExperiment"))
        x <- rowData(x)
    if (is(subject, "SummarizedExperiment"))
        subject <- rowData(subject)
    distanceToNearest(x, subject, ignore.strand=ignore.strand,
                      ...)
}

local({
    .signatures <- list(
        c("SummarizedExperiment", "ANY"),
        c("ANY", "SummarizedExperiment"),
        c("SummarizedExperiment", "SummarizedExperiment"))

    for (.sig in .signatures) {
        .funs <- c("nearest", "precede", "follow")
        tmpl <- function(x, subject, select = c("arbitrary", "all"),
                         ignore.strand = FALSE, ...) {}
        environment(tmpl) <- parent.frame(2)
        for (.fun in .funs) {
            body(tmpl) <- substitute({
                select <- match.arg(select)
                if (is(x, "SummarizedExperiment"))
                    x <- rowData(x)
                if (is(subject, "SummarizedExperiment"))
                    subject <- rowData(subject)
                FUN(x=x, subject=subject, select=select,
                    ignore.strand=ignore.strand, ...)
            }, list(FUN=as.symbol(.fun)))
            setMethod(.fun, .sig, tmpl)
        }
        setMethod("nearest", c("SummarizedExperiment", "missing"),
            .SummarizedExperiment.nearest.missing)
        setMethod("compare", .sig, .SummarizedExperiment.compare)
        setMethod("Compare", .sig, .SummarizedExperiment.Compare)
        setMethod("distance", .sig, .SummarizedExperiment.distance)
        setMethod("distanceToNearest", .sig,
            .SummarizedExperiment.distanceToNearest)
    }
})

## additional getters / setters

setReplaceMethod("strand", "SummarizedExperiment",
    function(x, ..., value)
{
    strand(rowData(x)) <- value
    x
})

setReplaceMethod("ranges", "SummarizedExperiment",
    function(x, ..., value)
{
    ranges(rowData(x)) <- value
    x
})

## order, rank, sort

setMethod("order", "SummarizedExperiment",
    function(..., na.last = TRUE, decreasing = FALSE)
{
    args <- lapply(list(...), rowData)
    do.call("order",
            c(args, list(na.last=na.last, decreasing=decreasing)))
})

setMethod("rank", "SummarizedExperiment",
    function (x, na.last = TRUE,
              ties.method = c("average", "first", "random", "max", "min"))
{
    ties.method <- match.arg(ties.method)
    rank(rowData(x), na.last=na.last, ties.method=ties.method)
})

setMethod("sort", "SummarizedExperiment",
    function(x, decreasing = FALSE, ...)
{
    x[order(rowData(x), decreasing=decreasing),]
})

## seqinfo (also seqlevels, genome, seqlevels<-, genome<-), seqinfo<-

setMethod(seqinfo, "SummarizedExperiment",
    function(x)
{
    seqinfo(rowData(x))
})

setReplaceMethod("seqinfo", "SummarizedExperiment",
    function (x, new2old = NULL, force = FALSE, value)
{
    if (!is(value, "Seqinfo")) 
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    dangling_seqlevels <-
        getDanglingSeqlevels(rowData(x), new2old = new2old,
                             force = force, seqlevels(value))
    if (length(dangling_seqlevels) != 0L) 
        x <- x[!(seqnames(x) %in% dangling_seqlevels)]
    rowData(x) <-
        update(rowData(x),
               seqnames = makeNewSeqnames(x, new2old, seqlevels(value)),
               seqinfo = value)
    if (is.character(msg <- .valid.SummarizedExperiment(x)))
        stop(msg)
    x
})
