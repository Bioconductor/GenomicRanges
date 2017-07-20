### =========================================================================
### GPos objects
### -------------------------------------------------------------------------
###


setClass("GPos",
    contains="GRanges",
    representation(
        ranges="IPos"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

setMethod("names", "GPos", function(x) NULL)

setReplaceMethod("names", "GPos",
    function(x, value)
    {
        if (!is.null(value))
            stop(class(x), " objects don't accept names")
        x
    }
)

setMethod("pos", "GPos", function(x) pos(ranges(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Collapse runs of "stitchable genomic ranges"
###
### 2 genomic ranges are "stitchable" if, in addition to be stitchable from
### an integer ranges point-of-view (see .stitch_Ranges() in
### IRanges/R/IPos-class.R for what that means), they are also on the same
### chromosome and strand.

### .stitch_GenomicRanges() below takes any GenomicRanges derivative and
### returns a GRanges object (so is NOT an endomorphism).
### Note that this transformation preserves 'sum(width(x))'.
### Also note that this is an "inter range transformation". However unlike
### range(), reduce(), gaps(), or disjoin(), its result depends on the order
### of the elements in the input vector. It's also idempotent like range(),
### reduce(), and disjoin() (gaps() is not).

### TODO: Define and export stitch() generic and method for Ranges objects
### in the IRanges package (in inter-range-methods.R). Then make
### .stitch_GenomicRanges() the "stitch" method for GenomicRanges objects and
### support the 'ignore.strand' argument.

.stitch_GenomicRanges <- function(x, drop.empty.ranges=FALSE)
{
    if (length(x) == 0L)
        return(granges(x))  # returning GRanges() would loose the seqinfo

    x_seqnames <- seqnames(x)
    x_strand <- strand(x)
    x_start <- start(x)
    x_end <- end(x)
    new_run <- x_seqnames[-1L] != x_seqnames[-length(x)] |
        x_strand[-1L] != x_strand[-length(x)] |
        Rle(x_start[-1L] != x_end[-length(x)] + 1L)
    new_run_idx <- which(new_run)
    start_idx <- c(1L, new_run_idx + 1L)
    end_idx <- c(new_run_idx, length(x))

    ans_ranges <- IRanges(x_start[start_idx], x_end[end_idx])

    if (drop.empty.ranges) {
        keep_idx <- which(width(ans_ranges) != 0L)
        ans_ranges <- ans_ranges[keep_idx]
        start_idx <- start_idx[keep_idx]
    }

    ans_seqnames <- x_seqnames[start_idx]
    ans_strand <- x_strand[start_idx]
    ans_mcols <- S4Vectors:::make_zero_col_DataFrame(length(start_idx))
    ans_seqinfo <- seqinfo(x)

    ## To be as fast as possible, we don't use internal low-level constructor
    ## new_GRanges() and we don't check the new object.
    new2("GRanges", seqnames=ans_seqnames,
                    ranges=ans_ranges,
                    strand=ans_strand,
                    elementMetadata=ans_mcols,
                    seqinfo=ans_seqinfo,
                    check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

### TODO


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### Note that if 'pos_runs' is a GPos instance with no metadata or metadata
### columns, then 'identical(GPos(pos_runs), pos_runs)' is TRUE.
GPos <- function(pos_runs=GRanges())
{
    if (!is(pos_runs, "GenomicRanges"))
        pos_runs <- as(pos_runs, "GenomicRanges", strict=FALSE)
    suppressWarnings(ans_len <- sum(width(pos_runs)))
    if (is.na(ans_len))
        stop("too many genomic positions in 'pos_runs'")
    ans_seqnames <- rep.int(seqnames(pos_runs), width(pos_runs))
    ans_ranges <- IPos(ranges(pos_runs))
    ans_strand <- rep.int(strand(pos_runs), width(pos_runs))
    ans_mcols <- S4Vectors:::make_zero_col_DataFrame(ans_len)
    new2("GPos", seqnames=ans_seqnames, ranges=ans_ranges, strand=ans_strand,
                 elementMetadata=ans_mcols, seqinfo=seqinfo(pos_runs),
                 check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

.from_GenomicRanges_to_GPos <- function(from)
{
    if (!all(width(from) == 1L))
        stop(wmsg("all the ranges in the ", class(from), " object to ",
                  "coerce to GPos must have a width of 1"))
    if (!is.null(names(from))) {
        names(from) <- NULL
        warning(wmsg("because a GPos object cannot hold them, the names ",
                     "on the ", class(from), " object couldn't be ",
                     "propagated during its coercion to GPos"))
    }
    from@ranges <- IPos(from@ranges)
    from
}
setAs("GenomicRanges", "GPos", .from_GenomicRanges_to_GPos)

### Automatic coercion method from GPos to GRanges silently returns
### a broken object (unfortunately these dummy automatic coercion methods
### don't bother to validate the object they return). So we overwrite it.
.from_GPos_to_GRanges <- function(from)
{
    class(from) <- "GRanges"  # temporarily broken GRanges instance!
    from@ranges <- as(from@ranges, "IRanges")  # now fixed :-)
    from
}
setAs("GPos", "GRanges", .from_GPos_to_GRanges)

### The "as.data.frame" method for GenomicRanges objects works on a GPos
### object but returns a data.frame with identical "start" and "end" columns,
### and a "width" column filled with 1. We overwrite it to return a data.frame
### with a "pos" column instead of the "start" and "end" columns, and no
### "width" column.
### TODO: Turn this into an S3/S4 combo for as.data.frame.GPos
setMethod("as.data.frame", "GPos",
    function(x, row.names=NULL, optional=FALSE, ...)
    {
        mcols_df <- as.data.frame(mcols(x), ...)
        data.frame(seqnames=as.factor(seqnames(x)),
                   pos=pos(x),
                   strand=as.factor(strand(x)),
                   mcols_df,
                   stringsAsFactors=FALSE)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

.make_naked_matrix_from_GPos <- function(x)
{
    x_len <- length(x)
    x_mcols <- mcols(x)
    x_nmc <- if (is.null(x_mcols)) 0L else ncol(x_mcols)
    ans <- cbind(seqnames=as.character(seqnames(x)),
                 pos=as.character(pos(x)),
                 strand=as.character(strand(x)))
    if (x_nmc > 0L) {
        tmp <- do.call(data.frame, c(lapply(x_mcols, showAsCell),
                                     list(check.names=FALSE)))
        ans <- cbind(ans, `|`=rep.int("|", x_len), as.matrix(tmp))
    }
    ans
}

show_GPos <- function(x, margin="",
                      print.classinfo=FALSE, print.seqinfo=FALSE)
{
    x_class <- class(x)
    x_len <- length(x)
    x_mcols <- mcols(x)
    x_nmc <- if (is.null(x_mcols)) 0L else ncol(x_mcols)
    cat(x_class, " object with ",
        x_len, " ", ifelse(x_len == 1L, "position", "positions"),
        " and ",
        x_nmc, " metadata ", ifelse(x_nmc == 1L, "column", "columns"),
        ":\n", sep="")
    ## S4Vectors:::makePrettyMatrixForCompactPrinting() assumes that head()
    ## and tail() work on 'xx'.
    xx <- as(x, "GPos")
    out <- S4Vectors:::makePrettyMatrixForCompactPrinting(xx,
                .make_naked_matrix_from_GPos)
    if (print.classinfo) {
        .COL2CLASS <- c(
            seqnames="Rle",
            pos="integer",
            strand="Rle"
        )
        classinfo <-
            S4Vectors:::makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
        ## A sanity check, but this should never happen!
        stopifnot(identical(colnames(classinfo), colnames(out)))
        out <- rbind(classinfo, out)
    }
    if (nrow(out) != 0L)
        rownames(out) <- paste0(margin, rownames(out))
    ## We set 'max' to 'length(out)' to avoid the getOption("max.print")
    ## limit that would typically be reached when 'showHeadLines' global
    ## option is set to Inf.
    print(out, quote=FALSE, right=TRUE, max=length(out))
    if (print.seqinfo) {
        cat(margin, "-------\n", sep="")
        cat(margin, "seqinfo: ", summary(seqinfo(x)), "\n", sep="")
    }
}

setMethod("show", "GPos",
    function(object)
        show_GPos(object, margin="  ",
                  print.classinfo=TRUE, print.seqinfo=TRUE)
)

