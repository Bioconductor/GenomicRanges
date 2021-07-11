### =========================================================================
### GPos objects
### -------------------------------------------------------------------------
###


setClass("GPos",
    contains=c("GenomicPos", "GRanges"),
    representation(
        "VIRTUAL",
        ranges="IPos"
    )
)

setClass("UnstitchedGPos",
    contains="GPos",
    representation(
        ranges="UnstitchedIPos"
    )
)

setClass("StitchedGPos",
    contains="GPos",
    representation(
        ranges="StitchedIPos"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

.OLD_GPOS_INSTANCE_MSG <- c(
    "Starting with BioC 3.10, the class attribute of all ",
    "GPos **instances** needs to be set to \"StitchedGPos\". ",
    "Please update this object with 'updateObject(object, verbose=TRUE)' ",
    "and re-serialize it."
)

.validate_GPos <- function(x)
{
    if (class(x) == "GPos")
        return(paste(.OLD_GPOS_INSTANCE_MSG, collapse=""))

    NULL
}

setValidity2("GPos", .validate_GPos)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("pos", "GPos", function(x) pos(ranges(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Collapse runs of "stitchable genomic ranges"
###
### 2 genomic ranges are "stitchable" if, in addition to be stitchable from
### an integer ranges point-of-view (see stitch_IntegerRanges() in
### IRanges/R/IPos-class.R for what that means), they are also on the same
### chromosome and strand.

### stitch_GenomicRanges() below takes any GenomicRanges derivative and
### returns a GRanges object (so is NOT an endomorphism). Note that this
### transformation preserves 'sum(width(x))'.
### Also note that this is an "inter range transformation". However unlike
### range(), reduce(), gaps(), or disjoin(), its result depends on the order
### of the elements in the input vector. It's also idempotent like range(),
### reduce(), and disjoin() (gaps() is not).

### TODO: Define and export stitch() generic and method for IntegerRanges
### objects in the IRanges package (in inter-range-methods.R). Then make
### stitch_GenomicRanges() and stitch_StitchedGPos() the "stitch" methods
### for GenomicRanges and StitchedGPos objects, respectively, and support
### the 'ignore.strand' argument.

### To be as fast as possible, we don't use internal low-level constructor
### new_GRanges() and we don't check the new object.
.new_stitched_GRanges <- function(seqnames, ranges, strand, seqinfo)
{
    mcols <- make_zero_col_DFrame(length(ranges))
    new2("GRanges", seqnames=seqnames,
                    ranges=ranges,
                    strand=strand,
                    elementMetadata=mcols,
                    seqinfo=seqinfo,
                    check=FALSE)
}

stitch_GenomicRanges <- function(x)
{
    if (length(x) == 0L)
        return(granges(x, use.names=FALSE))  # returning GRanges() would loose
                                             # the seqinfo

    x_seqnames <- seqnames(x)
    x_strand <- strand(x)
    x_start <- start(x)
    x_end <- end(x)

    ## Find runs of stitchable elements along 'x'.
    ## Each run is described by the indices of its first ('run_from') and
    ## last ('run_to') elements in 'x'.
    ## The runs form a partitioning of 'x'.
    is_new_run <- x_seqnames[-1L] != x_seqnames[-length(x)] |
                  x_strand[-1L] != x_strand[-length(x)] |
                  Rle(x_start[-1L] != x_end[-length(x)] + 1L)
    new_run_idx <- which(is_new_run)
    run_from <- c(1L, new_run_idx + 1L)
    run_to <- c(new_run_idx, length(x))

    ans_ranges <- IRanges(x_start[run_from], x_end[run_to])
    ans_seqnames <- x_seqnames[run_from]  # same as x_seqnames[run_to]
    ans_strand <- x_strand[run_from]      # same as x_strand[run_to]
    .new_stitched_GRanges(ans_seqnames, ans_ranges, ans_strand, seqinfo(x))
}

stitch_StitchedGPos <- function(x)
{
    if (length(x) == 0L)
        return(granges(x, use.names=FALSE))  # returning GRanges() would loose
                                             # the seqinfo

    x_seqnames <- seqnames(x)
    x_strand <- strand(x)

    ## Find runs of identical (seqnames, strand) pairs along 'x'.
    ## The runs are described by IRanges object 'runs'.
    ## They form a partitioning of 'x'.
    is_new_run <- x_seqnames[-1L] != x_seqnames[-length(x)] |
                  x_strand[-1L] != x_strand[-length(x)]
    new_run_idx <- which(is_new_run)
    run_from <- c(1L, new_run_idx + 1L)
    run_to <- c(new_run_idx, length(x))
    runs <- IRanges(run_from, run_to)

    ans_ranges <- IRanges:::extract_pos_runs_by_ranges(x@ranges@pos_runs, runs)
    breakpoints <- cumsum(width(ans_ranges))
    ans_seqnames <- x_seqnames[breakpoints]
    ans_strand <- x_strand[breakpoints]
    .new_stitched_GRanges(ans_seqnames, ans_ranges, ans_strand, seqinfo(x))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### High-level GPos constructor.
GPos <- function(seqnames=NULL, pos=NULL, strand=NULL,
                 ..., seqinfo=NULL, seqlengths=NULL, stitch=NA)
{
    mcols <- DataFrame(..., check.names=FALSE)

    if (!is.null(pos)) {
        pos <- IPos(pos, stitch=stitch)
    } else if (is.null(seqnames)) {
        pos <- IPos(stitch=stitch)
    } else {
        if (is(seqnames, "GPos")) {
            x <- seqnames
        } else {
            x <- as(seqnames, "GRanges")
        }
        x_ranges <- x@ranges  # either IPos or IRanges
        pos <- IPos(x_ranges, stitch=stitch)
        seqnames <- x@seqnames
        if (is(x_ranges, "IRanges"))  # i.e. 'x' is not a GPos
            seqnames <- rep.int(seqnames, width(x_ranges))
        if (is.null(strand)) {
            strand <- x@strand
            if (is(x_ranges, "IRanges"))  # i.e. 'x' is not a GPos
                strand <- rep.int(strand, width(x_ranges))
        }
        if (length(mcols) == 0L && is(x, "GPos"))
            mcols <- mcols(x, use.names=FALSE)
        if (is.null(seqinfo))
            seqinfo <- seqinfo(x)
    }

    seqinfo <- normarg_seqinfo2(seqinfo, seqlengths)

    ## We use as.character() to get rid of the "package" attribute on the
    ## class name returned by class(pos).
    Class <- sub("IPos$", "GPos", as.character(class(pos)))

    new_GRanges(Class, seqnames=seqnames, ranges=pos, strand=strand,
                       mcols=mcols, seqinfo=seqinfo)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

.try_to_coerce_to_GRanges_first <- function(from, to)
{
    if (is(from, "GRanges"))
        return(from)
    from <- try(as(from, "GRanges"), silent=TRUE)
    if (inherits(from, "try-error"))
        stop(wmsg("object to coerce to ", to, " ",
                  "couldn't be coerced to GRanges first"))
    from
}
.check_GenomicRanges_for_coercion_to_GPos <- function(from, to)
{
    if (!all(width(from) == 1L))
        stop(wmsg("all the ranges in the object to ",
                  "coerce to ", to, " must have a width of 1"))
}
.from_ANY_to_UnstitchedGPos <- function(from)
{
    from <- .try_to_coerce_to_GRanges_first(from, "UnstitchedGPos")
    .check_GenomicRanges_for_coercion_to_GPos(from, "UnstitchedGPos")
    class(from) <- "UnstitchedGPos"  # temporarily broken instance!
    from@ranges <- as(from@ranges, "UnstitchedIPos")  # now fixed :-)
    from
}
.from_ANY_to_StitchedGPos <- function(from)
{
    from <- .try_to_coerce_to_GRanges_first(from, "StitchedGPos")
    .check_GenomicRanges_for_coercion_to_GPos(from, "StitchedGPos")
    class(from) <- "StitchedGPos"  # temporarily broken instance!
    from@ranges <- as(from@ranges, "StitchedIPos")  # now fixed :-)
    from
}
setAs("ANY", "UnstitchedGPos", .from_ANY_to_UnstitchedGPos)
setAs("ANY", "StitchedGPos", .from_ANY_to_StitchedGPos)
setAs("ANY", "GPos", .from_ANY_to_UnstitchedGPos)

### Yes, we also need to define the 3 coercion methods below, even though
### they seem redundant with the 3 coercion methods above. This is because
### the oh-so-smart methods package wants to automatically define these
### coercion methods in case they are not explicitly defined by the user.
### Unfortunately, and not too surprisingly, these automatic coercion
### methods get it wrong! How could they possibly know what they are doing?
setAs("GRanges", "UnstitchedGPos", .from_ANY_to_UnstitchedGPos)
setAs("GRanges", "StitchedGPos", .from_ANY_to_StitchedGPos)
setAs("GRanges", "GPos", .from_ANY_to_UnstitchedGPos)

### Of course we want 'as(UnstitchedGPos, "GRanges", strict=FALSE)' and
### 'as(StitchedGPos, "GRanges", strict=FALSE)' to do the right thing
### (i.e. to be no-ops), but, unfortunately, as() won't do that if
### coerce,UnstitchedGPos,GRanges and coerce,StitchedGPos,GRanges methods
### are defined, because, in this case, as() will **always** call the
### method, EVEN WHEN strict=FALSE AND THE OBJECT TO COERCE ALREADY DERIVES
### FROM THE TARGET CLASS! (This is a serious flaw in as() current
### design/implementation but I wouldn't be surprised if someone argued
### that this is a feature and working as intended.)
### Anyway, a workaround is to support the 'strict=FALSE' case at the level
### of the coerce() method itself. However setAs() doesn't let us do that
### so this is why we use setMethod("coerce", ...) to define these methods.
from_GPos_to_GRanges <- function(from, to="GRanges", strict=TRUE)
{
    if (!isTRUEorFALSE(strict))
        stop("'strict' must be TRUE or FALSE")
    if (!strict)
        return(from)
    class(from) <- "GRanges"  # temporarily broken instance!
    from@ranges <- as(from@ranges, "IRanges")  # now fixed :-)
    from
}
setMethod("coerce", c("UnstitchedGPos", "GRanges"), from_GPos_to_GRanges)
setMethod("coerce", c("StitchedGPos", "GRanges"), from_GPos_to_GRanges)
### One might think that defining the coerce,GPos,GRanges method below would
### cover the UnstitchedGPos->GRanges and StitchedGPos->GRanges cases, but no
### such luck! Again, this is because the oh-so-smart methods package wants
### to automatically define the 2 coercion methods above in case they are not
### explicitly defined by the user. And once again, these automatic coercion
### methods get it wrong!
### For the same reason, UnstitchedGPos or StitchedGPos extensions (like the
### CTSS class in the CAGEr package) need to define a coercion method to
### GRanges otherwise they'll also get a broken automatic coercion method.
### They can do this with (from_GPos_to_GRanges is exported):
###   setMethod("coerce", c("CTSS", "GRanges"), from_GPos_to_GRanges)
#setMethod("coerce", c("GPos", "GRanges"), from_GPos_to_GRanges)

### S3/S4 combo for as.data.frame.GPos
### The "as.data.frame" method for GenomicRanges objects works on a GPos
### object but returns a data.frame with identical "start" and "end" columns,
### and a "width" column filled with 1. We overwrite it to return a data.frame
### with a "pos" column instead of the "start" and "end" columns, and no
### "width" column.
.as.data.frame.GPos <- function(x, row.names=NULL, optional=FALSE)
{
    if (!identical(optional, FALSE))
        warning(wmsg("'optional' argument was ignored"))
    x_mcols <- mcols(x, use.names=FALSE)  # always a DataFrame parallel to 'x'
    data.frame(seqnames=as.factor(seqnames(x)),
               pos=pos(x),
               strand=as.factor(strand(x)),
               as.data.frame(x_mcols),
               row.names=row.names,
               stringsAsFactors=FALSE)
}
as.data.frame.GPos <- function(x, row.names=NULL, optional=FALSE, ...)
    .as.data.frame.GPos(x, row.names=NULL, optional=FALSE, ...)
setMethod("as.data.frame", "GPos", .as.data.frame.GPos)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### updateObject()
###
### Internal representation of GPos objects has changed in GenomicRanges
### 1.29.10 (Bioc 3.6).
###

.get_GPos_version <- function(object)
{
    if (class(object) != "GPos")
        return("current")

    if (!.hasSlot(object, "pos_runs"))
        return(">= 1.29.10 and < 1.37.7")

    return("< 1.29.10")
}

setMethod("updateObject", "GPos",
    function(object, ..., verbose=FALSE)
    {
        version <- .get_GPos_version(object)
        if (.hasSlot(object, "pos_runs")) {
            if (verbose)
                message("[updateObject] ", class(object), " object ",
                        "uses internal representation from\n",
                        "[updateObject] GenomicRanges ", version, ". ",
                        "Updating it ... ", appendLF=FALSE)
            ans <- GPos(object@pos_runs)
            mcols(ans) <- mcols(object)
            metadata(ans) <- metadata(object)
            if (verbose)
                message("OK")
            return(ans)
        }
        if (class(object) == "GPos") {
            if (verbose)
                message("[updateObject] Settting class attribute of GPos ",
                        "instance to \"StitchedGPos\" ... ", appendLF=FALSE)
            class(object) <- class(new("StitchedGPos"))
            if (verbose)
                message("OK")
        }
        if (IRanges:::get_IPos_version(object@ranges) == "current" && verbose) {
            message("[updateObject] Internal representation of ",
                    class(object), " object is current.\n",
                    "[updateObject] Nothing to update.")
        }

        callNextMethod()
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Display
###

### S3/S4 combo for summary.GPos
summary.GPos <- summary.IPos
setMethod("summary", "GPos", summary.GPos)

.from_GPos_to_naked_character_matrix_for_display <- function(x)
{
    m <- cbind(seqnames=showAsCell(seqnames(x)),
               pos=showAsCell(pos(x)),
               strand=showAsCell(strand(x)))
    cbind_mcols_for_display(m, x)
}
setMethod("makeNakedCharacterMatrixForDisplay", "GPos",
    .from_GPos_to_naked_character_matrix_for_display
)

show_GPos <- function(x, margin="",
                      print.classinfo=FALSE, print.seqinfo=FALSE)
{
    version <- .get_GPos_version(x)
    if (version != "current")
        stop(c(wmsg("This ", class(x), " object uses internal representation ",
                    "from GenomicRanges ", version, ", and so needs to be ",
                    "updated before it can be displayed or used. ",
                    "Please update it with:"),
               "\n\n    object <- updateObject(object, verbose=TRUE)",
               "\n\n  and re-serialize it."))
    cat(margin, summary(x), ":\n", sep="")
    ## makePrettyMatrixForCompactPrinting() assumes that head() and tail()
    ## work on 'x'.
    out <- makePrettyMatrixForCompactPrinting(x)
    if (print.classinfo) {
        .COL2CLASS <- c(
            seqnames="Rle",
            pos="integer",
            strand="Rle"
        )
        classinfo <- makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
        ## A sanity check, but this should never happen!
        stopifnot(identical(colnames(classinfo), colnames(out)))
        out <- rbind(classinfo, out)
    }
    if (nrow(out) != 0L)
        rownames(out) <- paste0(margin, "  ", rownames(out))
    ## We set 'max' to 'length(out)' to avoid the getOption("max.print")
    ## limit that would typically be reached when 'showHeadLines' global
    ## option is set to Inf.
    print(out, quote=FALSE, right=TRUE, max=length(out))
    if (print.seqinfo) {
        cat(margin, "  -------\n", sep="")
        cat(margin, "  seqinfo: ", summary(seqinfo(x)), "\n", sep="")
    }
}

setMethod("show", "GPos",
    function(object)
        show_GPos(object, print.classinfo=TRUE, print.seqinfo=TRUE)
)

