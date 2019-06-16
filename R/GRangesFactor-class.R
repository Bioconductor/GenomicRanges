### =========================================================================
### GRangesFactor objects
### -------------------------------------------------------------------------
###
### Factor objects with GRanges levels.
###

setClass("GRangesFactor",
    contains="Factor",
    representation(
        levels="GRanges",
        elementMetadata="DataFrame"
    )
)

setMethod("FactorToClass", "GRanges", function(x) "GRangesFactor")


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###
### Factor() works but this specialized GRangesFactor constructor adds some
### minor convenience and extra checks.
###

GRangesFactor <- function(x, levels, index=NULL, ...)
{
    if (missing(x) && missing(levels) && missing(levels))
        return(Factor(levels=GRanges(), index=integer(0), ...))
    if (!(missing(x) || is(x, "GRanges")))
        stop(wmsg("'x' must be a GRanges object"))
    if (!(missing(levels) || is(levels, "GRanges")))
        stop(wmsg("'levels' must be a GRanges object"))
    Factor(x, levels, index=index, ...)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("seqnames", "GRangesFactor",
    function(x) extractROWS(seqnames(levels(x)), x@index)
)

setMethod("start", "GRangesFactor",
    function(x) extractROWS(start(levels(x)), x@index)
)

setMethod("end", "GRangesFactor",
    function(x) extractROWS(end(levels(x)), x@index)
)

setMethod("width", "GRangesFactor",
    function(x) extractROWS(width(levels(x)), x@index)
)

setMethod("pos", "GRangesFactor",
    function(x) extractROWS(pos(levels(x)), x@index)
)

setMethod("strand", "GRangesFactor",
    function(x, ...) extractROWS(strand(levels(x), ...), x@index)
)

setMethod("seqinfo", "GRangesFactor", function(x) seqinfo(levels(x)))

setMethod("granges", "GRangesFactor",
    function(x, use.names=TRUE, use.mcols=FALSE, ...)
        granges(unfactor(x), use.names=use.names, use.mcols=use.mcols, ...)
)

setMethod("ranges", "GRangesFactor",
    function(x, use.names=TRUE, use.mcols=FALSE, ...)
        ranges(unfactor(x), use.names=use.names, use.mcols=use.mcols, ...)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###
### In addition to the basic coercions defined for Factor objects, we
### support coercion back and forth between GRanges and GRangesFactor.
###

setAs("ANY", "GRangesFactor",
    function(from)
    {
        if (!is(from, "GRanges"))
            from <- as(from, "GRanges")
        as(from, "Factor")
    }
)

setAs("Factor", "GRanges",
    function(from)
    {
        ans <- unfactor(from)
        if (!is(ans, "GRanges"))
            ans <- as(ans, "GRanges")
        ans
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

.print_levels <- function(levels, margin="", max.print.levels=20L)
{
    nlevels <- length(levels)
    cat(margin, nlevels, " ", if (nlevels == 1L) "level" else "levels", sep="")
    if (nlevels == 0L) {
        cat("\n")
        return(invisible(NULL))
    }
    cat(":\n")
    if (nlevels > max.print.levels) {
        level_strings <- as.character(head(levels, n=max.print.levels-1L))
        level_strings <- c(unname(level_strings), "...")
    } else {
        level_strings <- unname(as.character(levels))
    }
    levels_margin <- paste0(margin, "  ")
    old_width <- getOption("width")
    new_width <- old_width - nchar(levels_margin)
    options(width=new_width)
    on.exit(options(width=old_width))
    out <- capture.output(print(level_strings, quote=FALSE))
    if (nlevels > max.print.levels) {
        footer <- paste0(" ... ", nlevels - max.print.levels + 1L,
                         " more levels ...")
        out <- c(out, footer)
    }
    options(width=old_width)
    out <- paste0(levels_margin, out)
    cat(out, sep="\n")
}

.show_GRangesFactor <- function(x, margin="",
                                print.classinfo=FALSE, print.seqinfo=FALSE)
{
    show_GenomicRanges(x, margin=margin,
                       print.classinfo=print.classinfo,
                       print.seqinfo=print.seqinfo)
    .print_levels(levels(x), margin=margin, 20L)
}

setMethod("show", "GRangesFactor",
    function(object)
        .show_GRangesFactor(object, print.classinfo=TRUE)
)

