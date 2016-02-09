### =========================================================================
### GPos objects
### -------------------------------------------------------------------------
###
### GPos is a container for storing a set of genomic *positions* i.e.
### genomic ranges of length 1. It's more memory-efficient than GRanges when
### the object contains long runs of adjacent positions.
###

setClass("GPos",
    contains="GenomicRanges",
    representation(
        pos_runs="GRanges",
        elementMetadata="DataFrame"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("length", "GPos", function(x) sum(width(x@pos_runs)))

setMethod("names", "GPos", function(x) NULL)

setReplaceMethod("names", "GPos",
    function(x, value)
    {
        if (!is.null(value))
            stop(class(x), " objects don't accept names")
        x
    }
)

setMethod("seqnames", "GPos",
    function(x) rep.int(seqnames(x@pos_runs), width(x@pos_runs))
)

setGeneric("pos", function(x) standardGeneric("pos"))
setMethod("pos", "GPos", function(x) as.integer(ranges(x@pos_runs)))
setMethod("start", "GPos", function(x) pos(x))
setMethod("end", "GPos", function(x) pos(x))
setMethod("width", "GPos", function(x) rep.int(1L, length(x)))
setMethod("ranges", "GPos", function(x) IRanges(pos(x), width=1L))

setMethod("strand", "GPos",
    function(x) rep.int(strand(x@pos_runs), width(x@pos_runs))
)

setMethod("seqinfo", "GPos", function(x) seqinfo(x@pos_runs))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### Merge adjacent ranges.
### Returns a GRanges object (NOT an endomorphism).
### Note that this transformation preserves 'sum(width(x))'.
.merge_adjacent_ranges <- function(x, drop.empty.ranges=FALSE)
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
    ans_mcols <- new("DataFrame", nrows=length(start_idx))
    ans_seqinfo <- seqinfo(x)

    ## To be as fast as possible, we don't use internal constructor
    ## newGRanges() and we don't check the new object.
    new2("GRanges", seqnames=ans_seqnames,
                    ranges=ans_ranges,
                    strand=ans_strand,
                    elementMetadata=ans_mcols,
                    seqinfo=ans_seqinfo,
                    check=FALSE)
}

### Note that if 'pos_runs' is a GPos instance with no metadata or metadata
### columns, then 'identical(GPos(pos_runs), pos_runs)' is TRUE.
GPos <- function(pos_runs=GRanges())
{
    if (!is(pos_runs, "GenomicRanges"))
        pos_runs <- as(pos_runs, "GenomicRanges", strict=FALSE)
    suppressWarnings(ans_len <- sum(width(pos_runs)))
    if (is.na(ans_len))
        stop("too many genomic positions in 'pos_runs'")
    ans_mcols <- new("DataFrame", nrows=ans_len)
    ans_pos_runs <- .merge_adjacent_ranges(pos_runs, drop.empty.ranges=TRUE)
    new2("GPos", pos_runs=ans_pos_runs,
                 elementMetadata=ans_mcols,
                 metadata=pos_runs@metadata,
                 check=FALSE)
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion
###

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
### Subsetting
###

setMethod("extractROWS", "GPos",
    function(x, i)
    {
        i <- normalizeSingleBracketSubscript(i, x, as.NSBS=TRUE)
        ## TODO: Maybe make this the coercion method from NSBS to Ranges.
        if (is(i, "RangesNSBS")) {
            ir <- i@subscript
            ir <- ir[width(ir) != 0L]
        } else {
            ir <- as(as.integer(i), "IRanges")
        }
        map <- S4Vectors:::map_ranges_to_runs(width(x@pos_runs),
                                              start(ir), width(ir))
        ## Because 'ir' has no zero-width ranges, 'spanned_nrun' cannot
        ## contain zeroes and so 'Ltrim' and 'Rtrim' cannot contain garbbage.
        offset_nrun <- map[[1L]]
        spanned_nrun <- map[[2L]]
        Ltrim <- map[[3L]]
        Rtrim <- map[[4L]]
        run_idx <- S4Vectors:::fancy_mseq(spanned_nrun, offset_nrun)
        new_pos_runs <- x@pos_runs[run_idx]
        if (length(run_idx) != 0L) {
            Rtrim_idx <- cumsum(spanned_nrun)
            Ltrim_idx <- c(1L, Rtrim_idx[-length(Rtrim_idx)] + 1L)
            trimmed_start <- start(new_pos_runs)[Ltrim_idx] + Ltrim
            trimmed_end <- end(new_pos_runs)[Rtrim_idx] - Rtrim
            start(new_pos_runs)[Ltrim_idx] <- trimmed_start
            end(new_pos_runs)[Rtrim_idx] <- trimmed_end
        }
        x@pos_runs <- .merge_adjacent_ranges(new_pos_runs)
        mcols(x) <- extractROWS(mcols(x), i)
        x
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining
###

### 'Class' must be "GPos" or the name of a concrete subclass of GPos.
### 'objects' must be a list of GPos objects.
### Returns an instance of class 'Class'.
combine_GPos_objects <- function(Class, objects,
                                 use.names=TRUE, ignore.mcols=FALSE)
{
    if (!isSingleString(Class))
        stop("'Class' must be a single character string")
    if (!extends(Class, "GPos"))
        stop("'Class' must be the name of a class that extends GPos")
    if (!is.list(objects))
        stop("'objects' must be a list")
    if (!isTRUEorFALSE(use.names))
        stop("'use.names' must be TRUE or FALSE")
    ### TODO: Support 'use.names=TRUE'.
    if (use.names)
        stop("'use.names=TRUE' is not supported yet")
    if (!isTRUEorFALSE(ignore.mcols))
        stop("'ignore.mcols' must be TRUE or FALSE")

    if (length(objects) != 0L) {
        ## TODO: Implement (in C) fast 'elementIsNull(objects)' in S4Vectors
        ## that does 'sapply(objects, is.null, USE.NAMES=FALSE)', and use it
        ## here.
        null_idx <- which(sapply(objects, is.null, USE.NAMES=FALSE))
        if (length(null_idx) != 0L)
            objects <- objects[-null_idx]
    }
    if (length(objects) == 0L)
        return(new(Class))

    ## TODO: Implement (in C) fast 'elementIs(objects, class)' in S4Vectors
    ## that does 'sapply(objects, is, class, USE.NAMES=FALSE)', and use it
    ## here. 'elementIs(objects, "NULL")' should work and be equivalent to
    ## 'elementIsNull(objects)'.
    if (!all(sapply(objects, is, Class, USE.NAMES=FALSE)))
        stop("the objects to combine must be ", Class, " objects (or NULLs)")
    objects_names <- names(objects)
    names(objects) <- NULL  # so lapply(objects, ...) below returns an
                            # unnamed list

    ## Combine "pos_runs" slots.
    pos_runs_slots <- lapply(objects, function(x) x@pos_runs)
    ## TODO: Use combine_GRanges_objects() here when it's available.
    ans_pos_runs <- do.call(c, pos_runs_slots)

    suppressWarnings(ans_len <- sum(width(ans_pos_runs)))
    if (is.na(ans_len))
        stop("too many genomic positions to combine")

    ## Combine "mcols" slots. We don't need to use fancy
    ## S4Vectors:::rbind_mcols() for this because the "mcols" slot of a
    ## GPos object is guaranteed to be a DataFrame.
    if (ignore.mcols) {
        ans_mcols <- new("DataFrame", nrows=ans_len)
    } else  {
        mcols_slots <- lapply(objects, function(x) x@elementMetadata)
        ## Will fail if not all the GPos objects in 'objects' have
        ## exactly the same metadata cols.
        ans_mcols <- do.call(rbind, mcols_slots)
    }

    ## Make 'ans' and return it.
    new2(Class, pos_runs=ans_pos_runs, elementMetadata=ans_mcols, check=FALSE)
}

setMethod("c", "GPos",
    function (x, ..., ignore.mcols=FALSE, recursive=FALSE)
    {
        if (!identical(recursive, FALSE))
            stop("\"c\" method for GPos objects ",
                 "does not support the 'recursive' argument")
        if (missing(x)) {
            objects <- list(...)
            x <- objects[[1L]]
        } else {
            objects <- list(x, ...)
        }
        combine_GPos_objects(class(x), objects,
                             use.names=FALSE,
                             ignore.mcols=ignore.mcols)
    }
)

